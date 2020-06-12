##' Weighted particle filter
##'
##' A sequential importance sampling (particle filter) algorithm.
##' Unlike in \code{pfilter}, resampling is performed only when triggered by
##' deficiency in the effective sample size.
##'
##' \bold{This function is experimental and should be considered in alpha stage.
##' Both interface and underlying algorithms may change without warning at any time.
##' Please explore the function and give feedback via the \pkg{pomp} Issues page.}
##'
##' @name wpfilter
##' @rdname wpfilter
##' @aliases wpfilter wpfilter,ANY-method wpfilter,missing-method
##' wpfilterd_pomp-class wpfilterd_pomp
##' @author Aaron A. King
##' @family elementary POMP methods
##' @family particle filter methods
##'
##' @include pomp_class.R pomp.R rprocess_spec.R dmeasure_spec.R pfilter.R
##' @importFrom stats setNames
##'
##' @inheritParams pfilter
##' @param trigger numeric; if the effective sample size becomes smaller than \code{trigger * Np}, resampling is triggered.
##' @param target numeric; target power.
##'
##' @return
##' An object of class \sQuote{wpfilterd_pomp}, which extends class \sQuote{pomp}.
##' Information can be extracted from this object using the methods documented below.
##' 
##' @section Methods:
##' \describe{
##' \item{\code{\link{logLik}}}{ the estimated log likelihood}
##' \item{\code{\link{cond.logLik}}}{ the estimated conditional log likelihood}
##' \item{\code{\link{eff.sample.size}}}{the (time-dependent) estimated effective sample size}
##' \item{\code{\link{as.data.frame}}}{ coerce to a data frame}
##' \item{\code{\link{plot}}}{diagnostic plots}
##' }
##' 
##' @references
##'
##' \Arulampalam2002
##'
NULL

setClass(
  "wpfilterd_pomp",
  contains="pomp",
  slots=c(
    trigger="numeric",
    target="numeric",
    eff.sample.size="numeric",
    cond.loglik="numeric",
    Np="integer",
    loglik="numeric"
  ),
  prototype=prototype(
    trigger=0.0,
    target=0.5,
    eff.sample.size=numeric(0),
    cond.loglik=numeric(0),
    Np=as.integer(NA),
    loglik=as.double(NA)
  )
)

setGeneric(
  "wpfilter",
  function (data, ...)
    standardGeneric("wpfilter")
)

setMethod(
  "wpfilter",
  signature=signature(data="missing"),
  definition=function (...) {
    reqd_arg("wpfilter","data")
  }
)

setMethod(
  "wpfilter",
  signature=signature(data="ANY"),
  definition=function (data, ...) {
    undef_method("wpfilter",data)
  }
)

##' @name wpfilter-data.frame
##' @aliases wpfilter,data.frame-method
##' @rdname wpfilter
##' @export
setMethod(
  "wpfilter",
  signature=signature(data="data.frame"),
  definition=function (
    data,
    Np,
    params, rinit, rprocess, dmeasure,
    trigger = 1, target = 0.5,
    ...,
    verbose = getOption("verbose", FALSE)) {

    tryCatch(
      wpfilter.internal(
        data,
        Np=Np,
        rinit=rinit,
        rprocess=rprocess,
        dmeasure=dmeasure,
        params=params,
        trigger=trigger,
        target=target,
        ...,
        verbose=verbose
      ),
      error = function (e) pStop("wpfilter",conditionMessage(e))
    )

  }
)

##' @name wpfilter-pomp
##' @aliases wpfilter,pomp-method
##' @rdname wpfilter
##' @export
setMethod(
  "wpfilter",
  signature=signature(data="pomp"),
  definition=function (
    data,
    Np,
    trigger = 1, target = 0.5,
    ...,
    verbose = getOption("verbose", FALSE)) {

    tryCatch(
      wpfilter.internal(
        data,
        Np=Np,
        trigger=trigger,
        target=target,
        ...,
        verbose=verbose
      ),
      error = function (e) pStop("wpfilter",conditionMessage(e))
    )

  }
)

##' @name wpfilter-wpfilterd_pomp
##' @aliases wpfilter,wpfilterd_pomp-method
##' @rdname wpfilter
##' @export
setMethod(
  "wpfilter",
  signature=signature(data="wpfilterd_pomp"),
  function (data, Np, trigger, target,
    ..., verbose = getOption("verbose", FALSE)) {

    if (missing(Np)) Np <- data@Np
    if (missing(trigger)) trigger <- data@trigger
    if (missing(target)) target <- data@target

    wpfilter(as(data,"pomp"),Np=Np,
      trigger=trigger,target=target,
      ...,verbose=verbose)

  }
)

wpfilter.internal <- function (object, Np, trigger, target, ...,
  .gnsi = TRUE, verbose = FALSE) {

  verbose <- as.logical(verbose)

  object <- pomp(object,...,verbose=verbose)

  if (undefined(object@rprocess) || undefined(object@dmeasure))
    pStop_(paste(sQuote(c("rprocess","dmeasure")),collapse=", ")," are needed basic components.")

  gnsi <- as.logical(.gnsi)

  params <- coef(object)
  times <- time(object,t0=TRUE)
  ntimes <- length(times)-1

  Np <- np_check(Np,ntimes)
    
  pompLoad(object,verbose=verbose)
  on.exit(pompUnload(object,verbose=verbose))

  trigger <- as.numeric(trigger)
  if (length(trigger)!=1L || !is.finite(trigger) || trigger < 0)
    pStop_(sQuote("trigger")," should be a non-negative scalar.")
  target <- as.numeric(target)
  if (length(target)!=1L || !is.finite(trigger) || target < 0 || target > 1)
    pStop_(sQuote("target")," should be a scalar in [0,1].")

  x <- rinit(object,params=params,nsim=Np[1L],.gnsi=gnsi)

  loglik <- rep(NA,ntimes)
  eff.sample.size <- numeric(ntimes)
  W <- numeric(Np[1L])

  for (nt in seq_len(ntimes)) { ## main loop

    ## advance the state variables according to the process model
    X <- rprocess(object,x0=x,t0=times[nt],times=times[nt+1],params=params,.gnsi=gnsi)

    ## density of Y_t | X_t
    w <- dmeasure(object,y=object@data[,nt,drop=FALSE],x=X,
      times=times[nt+1],params=params,log=TRUE,.gnsi=gnsi)
    gnsi <- FALSE

    ## compute effective sample size, log-likelihood.
    xx <- .Call(P_wpfilter_comps,X,params,W,w,trigger,target,Np[nt+1])

    ## the following is triggered by the first illegal weight value
    if (is.integer(xx)) {
      illegal_dmeasure_error(
        time=times[nt+1],
        loglik=w[xx],
        datvals=object@data[,nt],
        states=X[,xx,1L],
        params=params
      )
    }

    x <- xx$states
    params <- xx$params[,1L]
    W <- xx$weights
    loglik[nt] <- xx$loglik
    eff.sample.size[nt] <- xx$ess
    
  } ## end of main loop

  new(
    "wpfilterd_pomp",
    object,
    trigger=trigger,
    target=target,
    eff.sample.size=eff.sample.size,
    cond.loglik=loglik,
    Np=as.integer(Np),
    loglik=sum(loglik)
  )
}
