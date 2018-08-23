##' The Liu and West Bayesian particle filter
##'
##' Modified version of the Liu and West (2001) algorithm.
##'
##' \code{bsmc2} uses a version of the original algorithm, but discards the auxiliary particle filter.
##' The modification appears to give superior performance for the same amount of effort.
##'
##' Samples from the prior distribution are drawn using the \code{rprior} component.
##' This can depend on elements of \code{params}, i.e., some of the elements of \code{params} can be \dQuote{hyperparameters}.
##' \code{Np} draws are made from the prior distribution.
##'
##' @name bsmc2
##' @docType methods
##' @rdname bsmc2
##' @include pomp_class.R workhorses.R pomp.R
##' @aliases bsmc2 bsmc2,missing-method bsmc2,ANY-method
##' @family particle filter methods
##' @family \pkg{pomp} parameter estimation methods
##'
##' @importFrom mvtnorm rmvnorm
##' @importFrom stats median cov setNames
##'
##' @inheritParams pomp
##' @inheritParams pfilter
##'
##' @param Np number of particles
##'
##' @param smooth Kernel density smoothing parameter.
##' The compensating shrinkage factor will be \code{sqrt(1-smooth^2)}.
##' Thus, \code{smooth=0} means that no noise will be added to parameters.
##' The general recommendation is that the value of \code{smooth} should be chosen close to 0 (e.g., \code{shrink} ~ 0.1).
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
    pomp_stop("bsmc2",sQuote("data")," is a required argument")
  }
)

setMethod(
  "bsmc2",
  signature=signature(data="ANY"),
  definition=function (data, ...) {
    stop(sQuote("bsmc2")," is not defined for objects of class ",
      sQuote(class(data)),call.=FALSE)
  }
)

##' @name bsmc2-data.frame
##' @rdname bsmc2
##' @aliases bsmc2-data.frame bsmc2,data.frame-method
##' @export
setMethod(
  "bsmc2",
  signature=signature(data="data.frame"),
  definition = function (data,
    params, rprior, rinit, rprocess, rmeasure, partrans,
    Np, smooth = 0.1, tol = 1e-17, max.fail = 0, ...,
    verbose = getOption("verbose", FALSE)) {

    object <- tryCatch(
      pomp(data,rinit=rinit,rprocess=rprocess,rmeasure=rmeasure,
        rprior=rprior,partrans=partrans,params=params,...,verbose=verbose),
      error = function (e) pomp_stop("bsmc2",conditionMessage(e))
    )

    bsmc2(
      object,
      Np=Np,
      smooth=smooth,
      tol=tol,
      verbose=verbose,
      max.fail=max.fail
    )

  }
)

##' @name bsmc2-pomp
##' @rdname bsmc2
##' @aliases bsmc2-pomp bsmc2,pomp-method
##' @export
setMethod(
  "bsmc2",
  signature=signature(data="pomp"),
  definition = function (data, params, Np, smooth = 0.1, tol = 1e-17,
    max.fail = 0, ..., verbose = getOption("verbose", FALSE)) {

    tryCatch(
      bsmc2.internal(
        data,
        params=params,
        Np=Np,
        smooth=smooth,
        tol=tol,
        max.fail=max.fail,
        ...,
        verbose=verbose
      ),
      error = function (e) pomp_stop("bsmc2",conditionMessage(e))
    )

  }
)

##' @name plot-bsmcd_pomp
##' @aliases plot,bsmcd_pomp-method
##' @rdname plot
##'
##' @param thin integer; when the number of samples is very large, it can be helpful to plot a random subsample:
##' \code{thin} specifies the size of this subsample.
##'
##' @export
setMethod(
  "plot",
  signature=signature(x="bsmcd_pomp"),
  definition=function (x, pars, thin, ...) {
    if (missing(pars)) pars <- x@est
    pars <- as.character(pars)
    if (length(pars)<1) pomp_stop("plot","no parameters to plot.")
    if (missing(thin)) thin <- Inf
    bsmc.plot(
      prior=partrans(x,x@prior,dir="fromEst"),
      post=partrans(x,x@post,dir="fromEst"),
      pars=pars,
      thin=thin,
      ...
    )
  }
)

bsmc2.internal <- function (object, Np, smooth, tol, max.fail, ...,
  .getnativesymbolinfo = TRUE, verbose) {

  ep <- character(0)

  object <- tryCatch(
    pomp(object,...,verbose=verbose),
    error = function (e) pomp_stop(ep,conditionMessage(e))
  )

  gnsi <- as.logical(.getnativesymbolinfo)

  if (missing(Np))
    pomp_stop(ep,sQuote("Np")," must be specified.")
  if (!missing(Np) && (length(Np) > 1 || !is.finite(Np) || Np < 1))
    pomp_stop(ep,sQuote("Np")," must be a positive integer.")
  Np <- as.integer(Np)

  params <- parmat(coef(object),Np)

  smooth <- as.numeric(smooth)
  if ((length(smooth)!=1) || (!is.finite(smooth)) || (smooth>1) || (smooth<=0))
    pomp_stop(ep,sQuote("smooth")," must be a scalar in (0,1]")

  tol <- as.numeric(tol)
  if (length(tol) != 1 || !is.finite(tol) || tol < 0)
    pomp_stop(ep,sQuote("tol")," should be a small positive number.")

  hsq <- smooth^2             #  see Liu & West eq(10.3.12)
  shrink <- sqrt(1-hsq)       #  'a' parameter of Liu & West

  pompLoad(object,verbose=verbose)
  on.exit(pompUnload(object,verbose=verbose))

  params <- rprior(object,params=params,.getnativesymbolinfo=gnsi)

  ntimes <- length(time(object))
  npars <- nrow(params)
  paramnames <- rownames(params)
  prior <- params

  times <- time(object,t0=TRUE)
  x <- rinit(object,params=params,.getnativesymbolinfo=gnsi)
  nvars <- nrow(x)

  params <- partrans(object,params,dir="toEst",.getnativesymbolinfo=gnsi)

  estind <- which(apply(params,1,function(x) diff(range(x)) > 0))
  est <- paramnames[estind]
  nest <- length(est)
  if (nest < 1) pomp_stop("","no parameters to estimate")

  evidence <- as.numeric(rep(NA,ntimes))
  eff.sample.size <- as.numeric(rep(NA,ntimes))
  nfail <- 0L

  mu <- array(data=NA,dim=c(nvars,Np,1L))
  rownames(mu) <- rownames(x)
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
      rmvnorm(n=Np,mean=rep(0,nest),sigma=hsq*params.var,method="svd"),
      error = function (e)
        pomp_stop("rmvnorm",conditionMessage(e))
    )

    if (!all(is.finite(pert))) pomp_stop(ep,"extreme particle depletion")

    params[estind,] <- m[estind,]+t(pert)

    tparams <- partrans(object,params,dir="fromEst",.getnativesymbolinfo=gnsi)

    xpred <- rprocess(object,xstart=x,times=times[c(nt,nt+1)],params=tparams,
      offset=1L,.getnativesymbolinfo=gnsi)

    ## evaluate likelihood of observation given xpred (from L&W AGM (4))
    weights <- dmeasure(object,y=object@data[,nt,drop=FALSE],x=xpred,
      times=times[nt+1],params=tparams,.getnativesymbolinfo=gnsi)

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
        pomp_stop(ep,"too many filtering failures")
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

    if (verbose)
      cat("effective sample size =",round(eff.sample.size[nt],1),"\n")

    ## Matrix with samples (columns) from filtering distribution theta.t | Y.t
    if (!all.fail) {
      smp <- .Call(systematic_resampling,weights)
      x <- x[,smp,drop=FALSE]
      params[estind,] <- params[estind,smp,drop=FALSE]
    }

  }

  if (nfail>0)
    pomp_warn(ep,nfail," filtering ",
      ngettext(nfail,"failure","failures")," occurred.")

  ## replace parameters with point estimate (posterior median)
  coef(object,transform=TRUE) <- apply(params,1,median)

  new(
    "bsmcd_pomp",
    object,
    post=partrans(object,params,dir="fromEst",.getnativesymbolinfo=gnsi),
    prior=prior,
    est=as.character(est),
    eff.sample.size=eff.sample.size,
    smooth=smooth,
    nfail=as.integer(nfail),
    cond.log.evidence=evidence,
    log.evidence=sum(evidence)
  )
}

bsmc.plot <- function (prior, post, pars, thin, ...) {
  p1 <- sample.int(n=ncol(prior),size=min(thin,ncol(prior)))
  p2 <- sample.int(n=ncol(post),size=min(thin,ncol(post)))
  if (!all(pars %in% rownames(prior))) {
    missing <- which(!(pars%in%rownames(prior)))
    pomp_stop("plot","unrecognized parameters: ",paste(sQuote(pars[missing]),collapse=","))
  }
  prior <- t(prior[pars,,drop=FALSE])
  post <- t(post[pars,,drop=FALSE])
  all <- rbind(prior,post)

  scplot <- function (x, y, ...) { ## prior, posterior pairwise scatterplot
    op <- par(new=TRUE)
    on.exit(par(op))
    i <- which(x[1L]==all[1L,])
    j <- which(y[1L]==all[1L,])
    points(prior[p1,i],prior[p1,j],pch=20,col=rgb(0.85,0.85,0.85,0.1),
      xlim=range(all[,i]),ylim=range(all[,j]))
    points(post[p2,i],post[p2,j],pch=20,col=rgb(0,0,1,0.01))
  }

  dplot <- function (x, ...) { ## marginal posterior histogram
    i <- which(x[1L]==all[1L,])
    d1 <- density(prior[,i])
    d2 <- density(post[,i])
    usr <- par("usr")
    op <- par(usr=c(usr[c(1L,2L)],0,1.5*max(d1$y,d2$y)))
    on.exit(par(op))
    polygon(d1,col=rgb(0.85,0.85,0.85,0.5))
    polygon(d2,col=rgb(0,0,1,0.5))
  }

  if (length(pars) > 1) {
    pairs(all,labels=pars,panel=scplot,diag.panel=dplot)
  }  else {
    d1 <- density(prior[,1])
    d2 <- density(post[,1])
    usr <- par("usr")
    op <- par(usr=c(usr[c(1L,2L)],0,1.5*max(d1$y,d2$y)))
    on.exit(par(op))
    plot(range(all[,1]),range(c(0,d1$y,d2$y)),type='n',
      xlab=pars,ylab="density")
    polygon(d1,col=rgb(0.85,0.85,0.85,0.5))
    polygon(d2,col=rgb(0,0,1,0.5))
  }
}
