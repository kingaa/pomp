##' The Liu and West Bayesian particle filter
##'
##' Modified version of the Liu and West (2001) algorithm.
##'
##' @name bsmc2
##' @docType methods
##' @rdname bsmc2
##' @include pomp_class.R workhorses.R pomp.R
##' @aliases bsmc2 bsmc2,missing-method bsmc2,ANY-method
##  @family particle filter methods
##'
##' @details
##' There are two ways to specify the prior distribution of particles.
##' If \code{params} is unspecified or is a named vector, \code{Np} draws are made from the prior distribution, as specified by \code{\link{rprior}}.
##' Alternatively, \code{params} can be specified as an \code{npars} x \code{Np} matrix (with rownames).
##'
##' \code{bsmc2} uses a version of the original algorithm, but discards the auxiliary particle filter.
##' The modification appears to give superior performance for the same amount of effort.
##'
##' @inheritParams pomp
##' @inheritParams pfilter
##' @param params parameters
##' @param Np number of particles
##' @param est Names of the parameters that are to be estimated.  No updates will be made to the other parameters.
##' If \code{est} is not specified, all parameters for which there is variation in \code{params} will be estimated.
##' @param smooth Kernel density smoothing parameter.
##' The compensating shrinkage factor will be \code{sqrt(1-smooth^2)}.
##' Thus, \code{smooth=0} means that no noise will be added to parameters.
##' The general recommendation is that the value of \code{smooth} should be chosen close to 0 (e.g., \code{shrink} ~ 0.1).
##' @param transform logical;
##' if \code{TRUE}, the algorithm operates on the transformed scale.
##'
##' @return
##' An object of class \sQuote{bsmcd_pomp}.
##' The following methods are avaiable:
##' \describe{
##' \item{\code{\link[=plot,bsmcd_pomp-method]{plot}}}{produces diagnostic plots}
##' \item{as.data.frame}{puts the prior and posterior samples into a data frame}
##' }
##'
##' @author Michael Lavine, Matthew Ferrari, Aaron A. King, Edward L. Ionides
##'
##' @references
##' Liu, J. and M. West.
##' Combining Parameter and State Estimation in Simulation-Based Filtering.
##' In A. Doucet, N. de Freitas, and N. J. Gordon, editors,
##' Sequential Monte Carlo Methods in Practice, pages 197-224.
##' Springer, New York, 2001.
NULL

## In annotation, L&W AGM == Liu & West "A General Algorithm"
##
## params = the initial particles for the parameter values;
##          these should be drawn from the prior distribution for the parameters
## est = names of parameters to estimate; other parameters are not updated.
## smooth = parameter 'h' from AGM

setClass(
  "bsmcd_pomp",
  contains="pomp",
  slots=c(
    transform="logical",
    post="array",
    prior="array",
    est="character",
    eff.sample.size="numeric",
    smooth="numeric",
    nfail="integer",
    cond.log.evidence="numeric",
    log.evidence="numeric"
  )
)

setGeneric(
  "bsmc2",
  function (data, ...)
    standardGeneric("bsmc2")
)

setMethod(
  "bsmc2",
  signature=signature(data="missing"),
  definition=function (...) {
    stop("in ",sQuote("bsmc2"),": ",sQuote("data")," is a required argument",call.=FALSE)
  }
)

setMethod(
  "bsmc2",
  signature=signature(data="ANY"),
  definition=function (data, ...) {
    stop(sQuote("bsmc2")," is not defined when ",sQuote("data")," is of class ",sQuote(class(data)),call.=FALSE)
  }
)

##' @name bsmc2-data.frame
##' @rdname bsmc2
##' @aliases bsmc2-data.frame bsmc2,data.frame-method
setMethod(
  "bsmc2",
  signature=signature(data="data.frame"),
  definition = function (data, params, rprior, rinit, rprocess, rmeasure,
    Np, est, smooth = 0.1, tol = 1e-17, max.fail = 0, transform = FALSE, ...,
    verbose = getOption("verbose", FALSE)) {

    object <- pomp(data,rinit=rinit,rprocess=rprocess,
      rmeasure=rmeasure,dprior=dprior,...,verbose=verbose)

    bsmc2(
      object,
      params=params,
      Np=Np,
      est=est,
      smooth=smooth,
      tol=tol,
      verbose=verbose,
      max.fail=max.fail,
      transform=transform
    )

  }
)

##' @name bsmc2-pomp
##' @rdname bsmc2
##' @aliases bsmc2-pomp bsmc2,pomp-method
setMethod(
  "bsmc2",
  signature=signature(data="pomp"),
  definition = function (data, params, Np, est, smooth = 0.1, tol = 1e-17,
    max.fail = 0, transform = FALSE, ...,
    verbose = getOption("verbose", FALSE)) {

    bsmc2.internal(
      data,
      params=params,
      Np=Np,
      est=est,
      smooth=smooth,
      tol=tol,
      max.fail=max.fail,
      transform=transform,
      ...,
      verbose=verbose
    )

  }
)

bsmc2.internal <- function (object, params, Np, est, smooth, tol,
  max.fail, transform, .getnativesymbolinfo = TRUE,
  ..., verbose) {

  ep <- paste0("in ",sQuote("bsmc2"),": ")

  object <- pomp(object,...,verbose=verbose)

  gnsi <- as.logical(.getnativesymbolinfo)
  transform <- as.logical(transform)

  if (missing(params)) {
    params <- coef(object)
    if (missing(Np)) stop(ep,sQuote("Np")," is not specified.",call.=FALSE)
  } else if (is.matrix(params)) {
    if (Np != ncol(params))
      warning(ep,sQuote("Np")," is ignored when ",sQuote("params"),
        " is a matrix.",call.=FALSE)
    Np <- ncol(params)
  } else if (is.data.frame(params)) {
    params <- t(data.matrix(params))
    Np <- ncol(params)
  } else if (is.list(params)) {
    params <- unlist(params)
    if (missing(Np)) stop(ep,sQuote("Np")," is not specified.",call.=FALSE)
  } else if (is.null(params)) {
    params <- numeric(0)
    if (missing(Np)) stop(ep,sQuote("Np")," is not specified.",call.=FALSE)
  }

  Np <- as.integer(Np)

  if (length(Np) > 1 || !is.finite(Np) || Np < 1)
    stop(ep,sQuote("Np")," must be a positive integer.",call.=FALSE)

  if ((!is.matrix(params)) && (Np > 1)) {
    params <- tryCatch(
      rprior(object,params=parmat(params,Np),.getnativesymbolinfo=gnsi),
      error = function (e) {
        stop(ep,sQuote("rprior")," error: ",conditionMessage(e),call.=FALSE)
      }
    )
  }

  if (!is.numeric(params) || is.null(rownames(params)))
    stop(ep,sQuote("params")," should be suppled as a numeric matrix with rownames",
      " or a named numeric vector.",call.=FALSE)

  params <- as.matrix(params)

  if (transform)
    params <- partrans(object,params,dir="toEst",.getnativesymbolinfo=gnsi)

  ntimes <- length(time(object))
  npars <- nrow(params)
  paramnames <- rownames(params)
  prior <- params

  if (missing(est))
    est <- paramnames[apply(params,1,function(x)diff(range(x))>0)]
  estind <- match(est,paramnames)
  npars.est <- length(estind)
  if (any(is.na(estind))) {
    ind <- which(is.na(estind))
    stop(ep,"parameter(s) ",
      paste(sapply(est[ind],sQuote),collapse=","),
      " not found.",call.=FALSE)
  }

  if (npars.est<1)
    stop(ep,"no parameters to estimate",call.=FALSE)

  smooth <- as.numeric(smooth)
  if ((length(smooth)!=1) || (!is.finite(smooth)) || (smooth>1) || (smooth<=0))
    stop(ep,sQuote("smooth")," must be a scalar in (0,1]",call.=FALSE)

  hsq <- smooth^2             #  see Liu & West eq(10.3.12)
  shrink <- sqrt(1-hsq)       #  'a' parameter of Liu & West

  tol <- as.numeric(tol)
  if (length(tol) != 1 || !is.finite(tol) || tol < 0)
    stop(ep,sQuote("tol")," should be a small positive number.",call.=FALSE)

  pompLoad(object,verbose=verbose)
  on.exit(pompUnload(object,verbose=verbose))

  xstart <- rinit(
    object,
    params=if (transform) {
      partrans(object,params,dir="fromEst",.getnativesymbolinfo=gnsi)
    } else {
      params
    },
    .getnativesymbolinfo=gnsi
  )
  nvars <- nrow(xstart)

  times <- time(object,t0=TRUE)
  x <- xstart

  evidence <- as.numeric(rep(NA,ntimes))
  eff.sample.size <- as.numeric(rep(NA,ntimes))
  nfail <- 0L

  mu <- array(data=NA,dim=c(nvars,Np,1L))
  rownames(mu) <- rownames(xstart)
  m  <- array(data=NA,dim=c(npars,Np))
  rownames(m) <- rownames(params)

  for (nt in seq_len(ntimes)) {

    ## calculate particle means ; as per L&W AGM (1)
    params.mean <- apply(params,1L,mean)
    ## calculate particle covariances : as per L&W AGM (1)
    params.var  <- cov(t(params[estind,,drop=FALSE]))

    if (verbose) {
      cat("at step ",nt," (time = ",times[nt+1],")\n",sep="")
      print(
        rbind(
          prior.mean=params.mean[estind],
          prior.sd=sqrt(diag(params.var))
        )
      )
    }

    m <- shrink*params+(1-shrink)*params.mean

    ## sample new parameter vector as per L&W AGM (3) and Liu & West eq(3.2)
    pert <- tryCatch(
      rmvnorm(
        n=Np,
        mean=rep(0,npars.est),
        sigma=hsq*params.var,
        method="svd"
      ),
      error = function (e) {
        stop(ep,"in ",sQuote("rmvnorm"),": ",conditionMessage(e),call.=FALSE)
      }
    )
    if (!all(is.finite(pert)))
      stop(ep,"extreme particle depletion",call.=FALSE) # nocov

    params[estind,] <- m[estind,]+t(pert)

    if (transform)
      tparams <- partrans(object,params,dir="fromEst",.getnativesymbolinfo=gnsi)

    xpred <- rprocess(
      object,
      xstart=x,
      times=times[c(nt,nt+1)],
      params=if (transform) {
        tparams
      } else {
        params
      },
      offset=1L,
      .getnativesymbolinfo=gnsi
    )

    ## evaluate likelihood of observation given xpred (from L&W AGM (4))
    weights <- tryCatch(
      dmeasure(
        object,
        y=object@data[,nt,drop=FALSE],
        x=xpred,
        times=times[nt+1],
        params=if (transform) {
          tparams
        } else {
          params
        },
        .getnativesymbolinfo=gnsi
      ),
      error = function (e) {
        stop(ep,conditionMessage(e),call.=FALSE)
      }
    )

    gnsi <- FALSE  ## all native symbols have been looked up

    ## evaluate weights as per L&W AGM (5)

    storeForEvidence <- log(mean(weights))

    x[,] <- xpred

    ## test for failure to filter
    dim(weights) <- NULL   ### needed?
    failures <- ((weights<tol) | (!is.finite(weights))) # test for NA weights
    all.fail <- all(failures)
    if (all.fail) {                     # all particles are lost
      if (verbose) {
        message("filtering failure at time t = ",times[nt+1])
      }
      nfail <- nfail+1
      if (nfail > max.fail)
        stop(ep,"too many filtering failures",call.=FALSE)
      evidence[nt] <- log(tol)          # worst log-likelihood
      weights <- rep(1/Np,Np)
      eff.sample.size[nt] <- 0
    } else {                  # not all particles are lost
      ## compute log-likelihood
      evidence[nt] <- storeForEvidence
      weights[failures] <- 0
      weights <- weights/sum(weights)
      ## compute effective sample-size
      eff.sample.size[nt] <- 1/crossprod(weights)
    }

    if (verbose) {
      cat("effective sample size =",round(eff.sample.size[nt],1),"\n")
    }

    ## Matrix with samples (columns) from filtering distribution theta.t | Y.t
    if (!all.fail) {
      smp <- .Call(systematic_resampling,weights)
      x <- x[,smp,drop=FALSE]
      params[estind,] <- params[estind,smp,drop=FALSE]
    }

  }

  if (nfail>0)
    warning(
      ep,nfail,
      ngettext(
        nfail,
        msg1=" filtering failure occurred.",
        msg2=" filtering failures occurred."
      ),
      call.=FALSE
    )

  ## replace parameters with point estimate (posterior median)
  coef(object,transform=transform) <- apply(params,1,median)

  new(
    "bsmcd_pomp",
    object,
    transform=transform,
    post=params,
    prior=prior,
    est=as.character(est),
    eff.sample.size=eff.sample.size,
    smooth=smooth,
    nfail=as.integer(nfail),
    cond.log.evidence=evidence,
    log.evidence=sum(evidence)
  )
}
