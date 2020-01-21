##' Iterated filtering: maximum likelihood by iterated, perturbed Bayes maps
##'
##' An iterated filtering algorithm for estimating the parameters of a partially-observed Markov process.
##' Running \code{mif2} causes the algorithm to perform a specified number of particle-filter iterations.
##' At each iteration, the particle filter is performed on a perturbed version of the model, in which the parameters to be estimated are subjected to random perturbations at each observation.
##' This extra variability effectively smooths the likelihood surface and combats particle depletion by introducing diversity into particle population.
##' As the iterations progress, the magnitude of the perturbations is diminished according to a user-specified cooling schedule.
##' The algorithm is presented and justified in Ionides et al. (2015).
##'
##' @name mif2
##' @rdname mif2
##' @include pfilter.R workhorses.R pomp_class.R safecall.R continue.R
##' @aliases mif2 mif2,missing-method mif2,ANY-method
##' @author Aaron A. King, Edward L. Ionides, Dao Nguyen
##' @family particle filter methods
##' @family pomp parameter estimation methods
##'
##' @importFrom utils head
##' @importFrom stats  weighted.mean
##'
##' @inheritParams pomp
##' @inheritParams pfilter
##' @param Nmif The number of filtering iterations to perform.
##' @param Np the number of particles to use in filtering.
##' This may be specified as a single positive integer, in which case the same number of particles will be used at each timestep.
##' Alternatively, if one wishes the number of particles to vary across timestep, one may specify \code{Np} either as a vector of positive integers (of length \code{length(time(object))}) or as a function taking a positive integer argument.
##' In the latter case, \code{Np(n)} must be a single positive integer,
##' representing the number of particles to be used at the \code{n}-th timestep:
##' \code{Np(1)} is the number of particles to use going from \code{timezero(object)} to \code{time(object)[1]},
##' \code{Np(2)}, from \code{time(object)[1]} to \code{time(object)[2]},
##' and so on.
##' @param rw.sd specification of the magnitude of the random-walk perturbations that will be applied to some or all model parameters.
##' Parameters that are to be estimated should have positive perturbations specified here.
##' The specification is given using the \code{\link{rw.sd}} function, which creates a list of unevaluated expressions.
##' The latter are evaluated in a context where the model time variable is defined (as \code{time}).
##' The expression \code{ivp(s)} can be used in this context as shorthand for \preformatted{ifelse(time==time[1],s,0).}
##' Likewise, \code{ivp(s,lag)} is equivalent to \preformatted{ifelse(time==time[lag],s,0).}
##' See below for some examples.
##'
##' The perturbations that are applied are normally distributed with the specified s.d.
##' If parameter transformations have been supplied, then the perturbations are applied on the transformed (estimation) scale.
##' @param cooling.type,cooling.fraction.50 specifications for the cooling schedule,
##' i.e., the manner and rate with which the intensity of the parameter perturbations is reduced with successive filtering iterations.
##' \code{cooling.type} specifies the nature of the cooling schedule.
##' See below (under \dQuote{Specifying the perturbations}) for more detail.
##'
##' @return
##' Upon successful completion, \code{mif2} returns an object of class
##' \sQuote{mif2d_pomp}.
##'
##' @section Methods:
##' The following methods are available for such an object:
##' \describe{
##' \item{\code{\link{continue}}}{ picks up where \code{mif2} leaves off and performs more filtering iterations. }
##' \item{\code{\link{logLik}}}{ returns the so-called \dfn{mif log likelihood} which is the log likelihood of the perturbed model, not of the focal model itself.
##' To obtain the latter, it is advisable to run several \code{\link{pfilter}} operations on the result of a \code{mif2} computatation.}
##' \item{\code{\link{coef}}}{ extracts the point estimate }
##' \item{\code{\link{eff.sample.size}}}{ extracts the effective sample size of the final filtering iteration}
##' }
##' Various other methods can be applied, including all the methods applicable to a \code{\link[=pfilter]{pfilterd_pomp}} object and all other \pkg{pomp} estimation algorithms and diagnostic methods.
##'
##' @section Specifying the perturbations:
##' The \code{rw.sd} function simply returns a list containing its arguments as unevaluated expressions.
##' These are then evaluated in a context containing the model \code{time} variable.  This allows for easy specification of the structure of the perturbations that are to be applied.
##' For example,
##' \preformatted{
##'     rw.sd(a=0.05, b=ifelse(0.2,time==time[1],0),
##'           c=ivp(0.2), d=ifelse(time==time[13],0.2,0),
##'           e=ivp(0.2,lag=13), f=ifelse(time<23,0.02,0))
##' }
##' results in perturbations of parameter \code{a} with s.d. 0.05 at every time step, while parameters \code{b} and \code{c} both get perturbations of s.d. 0.2 only before the first observation.
##' Parameters \code{d} and \code{e}, by contrast, get perturbations of s.d.  0.2 only before the thirteenth observation.
##' Finally, parameter \code{f} gets a random perturbation of size 0.02 before every observation falling before \eqn{t=23}.
##'
##' On the \eqn{m}-th IF2 iteration, prior to time-point \eqn{n}, the \eqn{d}-th parameter is given a random increment normally distributed with mean \eqn{0} and standard deviation \eqn{c_{m,n} \sigma_{d,n}}{c[m,n] sigma[d,n]}, where \eqn{c} is the cooling schedule and \eqn{\sigma}{sigma} is specified using \code{rw.sd}, as described above.
##' Let \eqn{N} be the length of the time series and \eqn{\alpha=}{alpha=}\code{cooling.fraction.50}.
##' Then, when \code{cooling.type="geometric"}, we have \deqn{c_{m,n}=\alpha^{\frac{n-1+(m-1)N}{50N}}.}{c[m,n]=alpha^((n-1+(m-1)N)/(50N)).}
##' When \code{cooling.type="hyperbolic"}, we have \deqn{c_{m,n}=\frac{s+1}{s+n+(m-1)N},}{c[m,n]=(s+1)/(s+n+(m-1)N),} where \eqn{s} satisfies \deqn{\frac{s+1}{s+50N}=\alpha.}{(s+1)/(s+50N)=alpha.}
##' Thus, in either case, the perturbations at the end of 50 IF2 iterations are a fraction \eqn{\alpha}{alpha} smaller than they are at first.
##'
##' @inheritSection pfilter Filtering failures
##'
##' @references
##'
##' \Ionides2015
##' 
NULL

setClass(
  'mif2d_pomp',
  contains='pfilterd_pomp',
  slots=c(
    Nmif = 'integer',
    rw.sd = 'matrix',
    cooling.type = 'character',
    cooling.fraction.50 = 'numeric',
    traces = 'matrix'
  )
)

setGeneric(
  "mif2",
  function (data, ...)
    standardGeneric("mif2")
)

setMethod(
  "mif2",
  signature=signature(data="missing"),
  definition=function (...) {
    reqd_arg("mif2","data")
  }
)

setMethod(
  "mif2",
  signature=signature(data="ANY"),
  definition=function (data, ...) {
    undef_method("mif2",data)
  }
)

##' @name mif2-data.frame
##' @aliases mif2,data.frame-method
##' @rdname mif2
##' @export
setMethod(
  "mif2",
  signature=signature(data="data.frame"),
  definition = function (data,
    Nmif = 1, rw.sd,
    cooling.type = c("geometric", "hyperbolic"), cooling.fraction.50,
    Np, tol = 0, max.fail = Inf,
    params, rinit, rprocess, dmeasure, partrans,
    ..., verbose = getOption("verbose", FALSE)) {

    tryCatch(
      mif2.internal(
        data,
        Nmif=Nmif,
        rw.sd=rw.sd,
        cooling.type=match.arg(cooling.type),
        cooling.fraction.50=cooling.fraction.50,
        Np=Np,
        tol=tol,
        max.fail=max.fail,
        params=params,
        rinit=rinit,
        rprocess=rprocess,
        dmeasure=dmeasure,
        partrans=partrans,
        ...,
        verbose=verbose
      ),
      error = function (e) pStop("mif2",conditionMessage(e))
    )

  }
)

##' @name mif2-pomp
##' @aliases mif2,pomp-method
##' @rdname mif2
##' @export
setMethod(
  "mif2",
  signature=signature(data="pomp"),
  definition = function (data,
    Nmif = 1, rw.sd,
    cooling.type = c("geometric", "hyperbolic"), cooling.fraction.50,
    Np, tol = 0, max.fail = Inf,
    ..., verbose = getOption("verbose", FALSE)) {

    tryCatch(
      mif2.internal(
        data,
        Nmif=Nmif,
        rw.sd=rw.sd,
        cooling.type=match.arg(cooling.type),
        cooling.fraction.50=cooling.fraction.50,
        Np=Np,
        tol=tol,
        max.fail=max.fail,
        ...,
        verbose=verbose
      ),
      error = function (e) pStop("mif2",conditionMessage(e))
    )

  }
)

##' @name mif2-pfilterd_pomp
##' @aliases mif2,pfilterd_pomp-method
##' @rdname mif2
##' @export
setMethod(
  "mif2",
  signature=signature(data="pfilterd_pomp"),
  definition = function (data,
    Nmif = 1, Np, tol, max.fail = Inf,
    ..., verbose = getOption("verbose", FALSE)) {

    if (missing(Np)) Np <- data@Np
    if (missing(tol)) tol <- data@tol

    mif2(
      as(data,"pomp"),
      Nmif=Nmif,
      Np=Np,
      tol=tol,
      max.fail=max.fail,
      ...,
      verbose=verbose
    )
  }
)

##' @name mif2-mif2d_pomp
##' @aliases mif2,mif2d_pomp-method
##' @rdname mif2
##' @export
setMethod(
  "mif2",
  signature=signature(data="mif2d_pomp"),
  definition = function (data,
    Nmif, rw.sd,
    cooling.type, cooling.fraction.50,
    ..., verbose = getOption("verbose", FALSE)) {

    if (missing(Nmif)) Nmif <- data@Nmif
    if (missing(rw.sd)) rw.sd <- data@rw.sd
    if (missing(cooling.type)) cooling.type <- data@cooling.type
    if (missing(cooling.fraction.50)) cooling.fraction.50 <- data@cooling.fraction.50

    mif2(
      as(data,"pfilterd_pomp"),
      Nmif=Nmif,
      rw.sd=rw.sd,
      cooling.type=cooling.type,
      cooling.fraction.50=cooling.fraction.50,
      ...,
      verbose=verbose
    )

  }
)

##' @name continue-mif2d_pomp
##' @aliases continue,mif2d_pomp-method
##' @rdname continue
##'
##' @param Nmif positive integer; number of additional filtering iterations to perform
##'
##' @export
setMethod(
  "continue",
  signature=signature(object="mif2d_pomp"),
  definition = function (object, Nmif = 1, ...) {

    ndone <- object@Nmif

    obj <- mif2(object,Nmif=Nmif,...,
      .ndone=ndone,.paramMatrix=object@paramMatrix)

    object@traces[ndone+1,c('loglik','nfail')] <- obj@traces[1L,c('loglik','nfail')]
    obj@traces <- rbind(
      object@traces,
      obj@traces[-1L,colnames(object@traces)]
    )
    names(dimnames(obj@traces)) <- c("iteration","variable")
    obj@Nmif <- as.integer(ndone+Nmif)

    obj
  }
)

mif2.internal <- function (object, Nmif, rw.sd,
  cooling.type, cooling.fraction.50,
  Np, tol = 0, max.fail = Inf,
  ..., verbose,
  .ndone = 0L, .indices = integer(0), .paramMatrix = NULL,
  .gnsi = TRUE) {

  verbose <- as.logical(verbose)

  object <- pomp(object,...,verbose=verbose)

  if (undefined(object@rprocess) || undefined(object@dmeasure))
    pStop_(paste(sQuote(c("rprocess","dmeasure")),collapse=", ")," are needed basic components.")

  gnsi <- as.logical(.gnsi)

  if (length(Nmif) != 1 || !is.numeric(Nmif) || !is.finite(Nmif) || Nmif < 1)
    pStop_(sQuote("Nmif")," must be a positive integer.")
  Nmif <- as.integer(Nmif)

  if (is.null(.paramMatrix)) {
    start <- coef(object)
  } else {  ## if '.paramMatrix' is supplied, 'start' is ignored
    start <- apply(.paramMatrix,1L,mean)
  }

  ntimes <- length(time(object))

  if (is.null(Np)) {
    pStop_(sQuote("Np")," must be specified.")
  } else if (is.function(Np)) {
    Np <- tryCatch(
      vapply(seq_len(ntimes),Np,numeric(1)),
      error = function (e) {
        pStop_("if ",sQuote("Np"),
          " is a function, it must return a single positive integer.")
      }
    )
  } else if (!is.numeric(Np)) {
    pStop_(sQuote("Np"),
      " must be a number, a vector of numbers, or a function.")
  }

  if (length(Np) == 1) {
    Np <- rep(Np,times=ntimes)
  } else if (length(Np) > ntimes) {
    if (Np[1L] != Np[ntimes+1] || length(Np) > ntimes+1) {
      pWarn("mif2","Np[k] ignored for k > ",sQuote("length(time(object))"),".")
    }
    Np <- head(Np,ntimes)
  } else if (length(Np) < ntimes) {
    pStop_(sQuote("Np")," must have length 1 or ",
      sQuote("length(time(object))"),".")
  }

  if (!all(is.finite(Np)) || any(Np <= 0))
    pStop_(sQuote("Np")," must be a positive integer.")

  Np <- as.integer(Np)
  Np <- c(Np,Np[1L])

  if (missing(rw.sd))
    pStop_(sQuote("rw.sd")," must be specified!")
  rw.sd <- perturbn.kernel.sd(rw.sd,time=time(object),paramnames=names(start))

  if (missing(cooling.fraction.50))
    pStop_(sQuote("cooling.fraction.50")," is a required argument.")
  if (length(cooling.fraction.50) != 1 || !is.numeric(cooling.fraction.50) ||
      !is.finite(cooling.fraction.50) || cooling.fraction.50 <= 0 ||
      cooling.fraction.50 > 1)
    pStop_(sQuote("cooling.fraction.50")," must be in (0,1].")
  cooling.fraction.50 <- as.numeric(cooling.fraction.50)

  cooling.fn <- mif2.cooling(
    type=cooling.type,
    fraction=cooling.fraction.50,
    ntimes=length(time(object))
  )

  if (is.null(.paramMatrix)) {
    paramMatrix <- array(data=start,dim=c(length(start),Np[1L]),
      dimnames=list(variable=names(start),rep=NULL))
  } else {
    paramMatrix <- .paramMatrix
  }

  traces <- array(dim=c(Nmif+1,length(start)+2),
    dimnames=list(iteration=seq.int(.ndone,.ndone+Nmif),
      variable=c('loglik','nfail',names(start))))
  traces[1L,] <- c(NA,NA,start)

  pompLoad(object,verbose=verbose)
  on.exit(pompUnload(object,verbose=verbose))

  paramMatrix <- partrans(object,paramMatrix,dir="toEst",
    .gnsi=gnsi)

  ## iterate the filtering
  for (n in seq_len(Nmif)) {

    pfp <- mif2.pfilter(
      object=object,
      params=paramMatrix,
      Np=Np,
      mifiter=.ndone+n,
      cooling.fn=cooling.fn,
      rw.sd=rw.sd,
      tol=tol,
      max.fail=max.fail,
      verbose=verbose,
      .indices=.indices,
      .gnsi=gnsi
    )

    gnsi <- FALSE

    paramMatrix <- pfp@paramMatrix
    traces[n+1,-c(1,2)] <- coef(pfp)
    traces[n,c(1,2)] <- c(pfp@loglik,pfp@nfail)
    .indices <- pfp@indices

    if (verbose) cat("mif2 iteration",n,"of",Nmif,"completed\n")

  }

  pfp@paramMatrix <- partrans(object,paramMatrix,dir="fromEst",
    .gnsi=gnsi)

  new(
    "mif2d_pomp",
    pfp,
    Nmif=Nmif,
    rw.sd=rw.sd,
    cooling.type=cooling.type,
    cooling.fraction.50=cooling.fraction.50,
    traces=traces
  )

}

mif2.cooling <- function (type, fraction, ntimes) {
  switch(
    type,
    geometric={
      factor <- fraction^(1/50)
      function (nt, m) {
        alpha <- factor^(nt/ntimes+m-1)
        list(alpha=alpha,gamma=alpha^2)
      }
    },
    hyperbolic={
      if (fraction < 1) {
        scal <- (50*ntimes*fraction-1)/(1-fraction)
        function (nt, m) {
          alpha <- (1+scal)/(scal+nt+ntimes*(m-1))
          list(alpha=alpha,gamma=alpha^2)
        }
      } else {
        function (nt, m) {
          list(alpha=1,gamma=1)
        }
      }
    }
  )
}

mif2.pfilter <- function (object, params, Np, mifiter, rw.sd, cooling.fn,
  tol = 0, max.fail = Inf, verbose, .indices = integer(0),
  .gnsi = TRUE) {

  tol <- as.numeric(tol)
  gnsi <- as.logical(.gnsi)
  verbose <- as.logical(verbose)
  mifiter <- as.integer(mifiter)
  Np <- as.integer(Np)

  if (length(tol) != 1 || !is.finite(tol) || tol < 0)
    pStop_(sQuote("tol")," should be a small nonnegative number.")

  if (tol != 0) {
    pWarn("mif2","the ",sQuote("tol")," argument is deprecated and will be removed in a future release. ","In the current release, the default value of ",
      sQuote("tol")," is 0. In future releases, the option to choose otherwise will be removed.")
  }

  do_ta <- length(.indices)>0L
  if (do_ta && length(.indices)!=Np[1L])
    pStop_(sQuote(".indices")," has improper length.")

  times <- time(object,t0=TRUE)
  ntimes <- length(times)-1

  loglik <- rep(NA,ntimes)
  eff.sample.size <- numeric(ntimes)
  nfail <- 0

  for (nt in seq_len(ntimes)) {

    ## perturb parameters
    pmag <- cooling.fn(nt,mifiter)$alpha*rw.sd[,nt]
    params <- .Call(P_randwalk_perturbation,params,pmag)

    tparams <- partrans(object,params,dir="fromEst",.gnsi=gnsi)

    ## get initial states
    if (nt == 1L) {
      x <- rinit(object,params=tparams)
    }

    ## advance the state variables according to the process model
    X <- rprocess(
      object,
      x0=x,
      t0=times[nt],
      times=times[nt+1],
      params=tparams,
      .gnsi=gnsi
    )

    ## determine the weights
    weights <-dmeasure(
      object,
      y=object@data[,nt,drop=FALSE],
      x=X,
      times=times[nt+1],
      params=tparams,
      log=FALSE,
      .gnsi=gnsi
    )

    if (!all(is.finite(weights))) {
      first <- which(!is.finite(weights))[1L]
      datvals <- object@data[,nt]
      weight <- weights[first]
      states <- X[,first,1L]
      params <- tparams[,first]
      msg <- nonfinite_dmeasure_error(time=times[nt+1],lik=weight,datvals,
        states,params)
      pStop_(msg)
    }
    gnsi <- FALSE

    ## compute weighted mean at last timestep
    if (nt == ntimes) {
      if (any(weights>0)) {
        coef(object,transform=TRUE) <- apply(params,1L,weighted.mean,w=weights)
      } else {
        pWarn("mif2","filtering failure at last filter iteration; using ",
          "unweighted mean for point estimate.")
        coef(object,transform=TRUE) <- apply(params,1L,mean)
      }
    }

    ## compute effective sample size, log-likelihood
    ## also do resampling if filtering has not failed
    xx <- .Call(P_pfilter_computations,x=X,params=params,Np=Np[nt+1],
      predmean=FALSE,predvar=FALSE,filtmean=FALSE,trackancestry=do_ta,
      doparRS=TRUE,weights=weights,tol=tol)

    all.fail <- xx$fail
    loglik[nt] <- xx$loglik
    eff.sample.size[nt] <- xx$ess
    if (do_ta) .indices <- .indices[xx$ancestry]

    x <- xx$states
    params <- xx$params

    if (all.fail) { ## all particles are lost
      nfail <- nfail+1
      if (verbose)
        message("filtering failure at time t = ",times[nt+1],".")
      if (nfail>max.fail)
        pStop_("too many filtering failures.")
    }

    if (verbose && (nt%%5==0))
      cat("mif2 pfilter timestep",nt,"of",ntimes,"finished.\n")

  }

  if (nfail>0)
    pWarn("mif2",nfail," filtering failure",ngettext(nfail,msg1="",msg2="s"),
      " occurred.")

  new(
    "pfilterd_pomp",
    as(object,"pomp"),
    paramMatrix=params,
    eff.sample.size=eff.sample.size,
    cond.loglik=loglik,
    indices=.indices,
    Np=Np,
    tol=tol,
    nfail=as.integer(nfail),
    loglik=sum(loglik)
  )
}

perturbn.kernel.sd <- function (rw.sd, time, paramnames) {

  if (is.matrix(rw.sd)) return(rw.sd)
  if (is(rw.sd,"safecall")) {
    enclos <- rw.sd@envir
    rw.sd <- as.list(rw.sd@call)[-1L]
  } else {
    pStop_(sQuote("rw.sd")," should be specified using the ",sQuote("rw.sd"),
      " function. See ",sQuote("?mif2"),".")
  }
  if (is.null(names(rw.sd)) | any(names(rw.sd)==""))
    pStop("rw.sd","parameters must be referenced by name.")
  if (!all(names(rw.sd) %in% paramnames)) {
    unrec <- names(rw.sd)[!names(rw.sd) %in% paramnames]
    pStop_("the following parameter(s), ",
      "given random walks in ",sQuote("rw.sd"),", are not present in ",
      sQuote("params"),": ",paste(sapply(unrec,sQuote),collapse=","),".")
  }
  ivp <- function (sd, lag = 1L) {
    sd*(seq_along(time)==lag)
  }
  sds <- lapply(rw.sd,eval,envir=list(time=time,ivp=ivp),enclos=enclos)
  for (n in names(sds)) {
    len <- length(sds[[n]])
    if (len==1) {
      sds[[n]] <- rep(sds[[n]],length(time))
    } else if (len!=length(time)) {
      pStop_(sQuote("rw.sd")," spec for parameter ",sQuote(n),
        " does not evaluate to a vector of the correct length (",
        sQuote("length(time(object))"),"=",length(time),").")
    }
  }
  do.call(rbind,sds)
}

##' rw.sd
##'
##' Specifying random-walk intensities.
##'
##' See \code{\link{mif2}} for details.
##'
##' @name rw.sd
##' @rdname rw_sd
##' @param \dots Specification of the random-walk intensities (as standard deviations).
##' @seealso \code{\link{mif2}}
##'
##' @export
rw.sd <- safecall
