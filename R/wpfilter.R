##' Weighted particle filter
##'
##' A weighted sequential Monte Carlo (particle filter) algorithm.
##' Resampling is performed according to a trigger.
##'
##' THIS FUNCTION IS EXPERIMENTAL.
##' IT MAY CHANGE WITHOUT WARNING AT ANY TIME.
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
##'
##' @return
##' An object of class \sQuote{wpfilterd_pomp}, which extends class \sQuote{pomp}.
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
    eff.sample.size="numeric",
    cond.loglik="numeric",
    Np="integer",
    loglik="numeric"
  ),
  prototype=prototype(
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
    ...,
    verbose = getOption("verbose", FALSE)) {

    tryCatch(
      wpfilter.internal(
        data,
        Np=Np,
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
  function (data, Np,
    ..., verbose = getOption("verbose", FALSE)) {

    if (missing(Np)) Np <- data@Np

    wpfilter(as(data,"pomp"),Np=Np,...,verbose=verbose)

  }
)

wpfilter.internal <- function (object, Np, ...,
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

  x <- rinit(object,params=params,nsim=Np[1L],.gnsi=gnsi)

  loglik <- rep(NA,ntimes)
  eff.sample.size <- numeric(ntimes)

  for (nt in seq_len(ntimes)) { ## main loop

    ## advance the state variables according to the process model
    X <- rprocess(object,x0=x,t0=times[nt],times=times[nt+1],params=params,.gnsi=gnsi)

    ## determine the weights
    weights <- dmeasure(object,y=object@data[,nt,drop=FALSE],x=X,
      times=times[nt+1],params=params,log=TRUE,.gnsi=gnsi)
    gnsi <- FALSE

    ## compute effective sample size, log-likelihood.
    ## also do resampling.
    xx <- .Call(P_wpfilter_comps,x=X,params=params,Np=Np[nt+1],weights=weights)

    ## the following is triggered by the first illegal weight value
    if (is.integer(xx)) {
      illegal_dmeasure_error(
        time=times[nt+1],
        loglik=weights[xx],
        datvals=object@data[,nt],
        states=X[,xx,1L],
        params=params
      )
    }

    loglik[nt] <- xx$loglik
    eff.sample.size[nt] <- xx$ess

    x <- xx$states
    params <- xx$params[,1L]

  } ## end of main loop

  new(
    "wpfilterd_pomp",
    object,
    eff.sample.size=eff.sample.size,
    cond.loglik=loglik,
    Np=as.integer(Np),
    loglik=sum(loglik)
  )
}
