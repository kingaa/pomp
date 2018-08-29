##' Spectrum matching
##'
##' Estimation of parameters by matching power spectra
##'
##' In spectrum matching, one attempts to minimize the discrepancy between a \acronym{POMP} model's predictions and data, as measured in the frequency domain by the power spectrum.
##'
##' \code{spect.match.objfun} constructs an objective function that measures the discrepancy.
##' It can be passed to any one of a variety of numerical optimization routines, which will adjust model parameters to minimize the discrepancies between the power spectrum of model simulations and that of the data.
##'
##' @name spect.match
##' @docType methods
##' @rdname spect_match
##' @family \pkg{pomp} parameter estimation methods
##' @aliases spect.match.objfun spect.match.objfun,missing-method spect.match.objfun,ANY-method
##'
##' @include spect.R probe_match.R loglik.R
#'
##' @inheritParams pomp
##' @inheritParams probe.match
##' @inheritParams spect
##'
##' @param weights optional numeric or function.
##' The mismatch between model and data is measured by a weighted average of mismatch at each frequency.
##' By default, all frequencies are weighted equally.
##' \code{weights} can be specified either as a vector (which must have length equal to the number of frequencies) or as a function of frequency.
##' If the latter, \code{weights(freq)} must return a nonnegative weight for each frequency.
##'
##' @return
##' \code{spect.match.objfun} constructs a stateful objective function for spectrum matching.
##' Specifically, \code{spect.match.objfun} returns an object of class \sQuote{spect_match_objfun}, which is a function suitable for use in an \code{\link{optim}}-like optimizer.
##' This function takes a single numeric-vector argument that is assumed to contain the parameters named in \code{est}, in that order.
##' When called, it will return the (optionally weighted) \eqn{L^2}{L2} distance between the data spectrum and simulated spectra.
##' It is a stateful function:
##' Each time it is called, it will remember the values of the parameters and the discrepancy measure.
##'
##' @inheritSection objfun Important Note
##'
##' @seealso \code{\link{spect}} \code{\link{optim}}
##' \code{\link[subplex]{subplex}} \code{\link[nloptr]{nloptr}}
##'
NULL

setClass(
  "spect_match_objfun",
  contains="function",
  slots=c(
    env="environment"
  )
)

setGeneric(
  "spect.match.objfun",
  function (data, ...)
    standardGeneric("spect.match.objfun")
)

setMethod(
  "spect.match.objfun",
  signature=signature(data="missing"),
  definition=function (...) {
    reqd_arg("spect.match.objfun","data")
  }
)

setMethod(
  "spect.match.objfun",
  signature=signature(data="ANY"),
  definition=function (data, ...) {
    undef_method("spect.match.objfun",data)
  }
)

##' @name spect.match.objfun-data.frame
##' @aliases spect.match.objfun,data.frame-method
##' @rdname spect_match
##' @export
setMethod(
  "spect.match.objfun",
  signature=signature(data="data.frame"),
  definition=function(data,
    rinit, rprocess, rmeasure, partrans, params,
    est = character(0), vars, nsim, seed = NULL, kernel.width,
    transform.data = identity,
    detrend = c("none","mean","linear","quadratic"),
    weights = 1, fail.value = NA, ...,
    verbose = getOption("verbose", FALSE)) {

    data <- tryCatch(
      pomp(data,rinit=rinit,rprocess=rprocess,rmeasure=rmeasure,
        partrans=partrans,params=params,...,verbose=verbose),
      error = function (e)
        pStop("spect.match.objfun",conditionMessage(e))
    )

    spect.match.objfun(data,est=est,vars=vars,nsim=nsim,seed=seed,
      kernel.width=kernel.width,transform.data=transform.data,detrend=detrend,
      weights=weights,fail.value=fail.value,verbose=verbose)

  }
)

##' @name spect.match.objfun-pomp
##' @aliases spect.match.objfun,pomp-method
##' @rdname spect_match
##' @export
setMethod(
  "spect.match.objfun",
  signature=signature(data="pomp"),
  definition=function(data, est = character(0),
    vars, nsim, seed = NULL, kernel.width, transform.data,
    detrend, weights = 1, fail.value = NA, ...,
    verbose = getOption("verbose", FALSE)) {

    smof.internal(data,est=est,vars=vars,nsim=nsim,seed=seed,
      kernel.width=kernel.width,transform.data=transform.data,detrend=detrend,
      weights=weights,fail.value=fail.value,...,verbose=verbose)

  }
)

##' @name spect.match.objfun-spectd_pomp
##' @aliases spect.match.objfun,spectd_pomp-method
##' @rdname spect_match
##' @export
setMethod(
  "spect.match.objfun",
  signature=signature(data="spectd_pomp"),
  definition=function(data, est = character(0), vars, nsim, seed = NULL,
    kernel.width, transform.data, detrend,
    weights = 1, fail.value = NA, ...,
    verbose = getOption("verbose", FALSE)) {

    if (missing(vars)) vars <- data@vars
    if (missing(nsim)) nsim <- data@nsim
    if (missing(kernel.width)) kernel.width <- data@kernel.width
    if (missing(transform.data)) transform.data <- data@transform.data
    if (missing(detrend)) detrend <- data@detrend

    spect.match.objfun(as(data,"pomp"),est=est,vars=vars,nsim=nsim,seed=seed,
      kernel.width=kernel.width,transform.data=transform.data,detrend=detrend,
      weights=weights,fail.value=fail.value,...,verbose=verbose)

  }
)

##' @name spect.match.objfun-spect_match_objfun
##' @aliases spect.match.objfun,spect_match_objfun-method
##' @rdname spect_match
##' @export
setMethod(
  "spect.match.objfun",
  signature=signature(data="spect_match_objfun"),
  definition=function(data, ..., verbose = getOption("verbose", FALSE)) {

    spect.match.objfun(data@env$object,...,verbose=verbose)

  }
)

smof.internal <- function (object, est, vars, nsim, seed,
  kernel.width, transform.data, detrend, weights, fail.value, ...) {

  object <- spect(object,vars=vars,nsim=nsim,seed=seed,...,
    kernel.width=kernel.width,transform.data=transform.data,detrend=detrend)

  fail.value <- as.numeric(fail.value)

  if (is.numeric(weights)) {
    if (length(weights)==1) {
      weights <- rep(weights,length(object@freq))
    } else if ((length(weights) != length(object@freq)))
      pStop_("if ",sQuote("weights"),
        " is provided as a vector, it must have length ",
        length(object@freq))
  } else if (is.function(weights)) {
    weights <- tryCatch(
      vapply(object@freq,weights,numeric(1)),
      error = function (e) {
        pStop_("problem with ",sQuote("weights")," function: ",
          conditionMessage(e))
      }
    )
  } else {
    pStop_(sQuote("weights"),
      " must be specified as a vector or as a function")
  }

  if (any(!is.finite(weights) | weights<0))
    pStop_(sQuote("weights")," should be nonnegative and finite")
  weights <- weights/mean(weights)


  params <- coef(object,transform=TRUE)

  if (missing(est)) est <- character(0)
  est <- as.character(est)
  est <- est[nzchar(est)]

  idx <- match(est,names(params))
  if (any(is.na(idx))) {
    missing <- est[is.na(idx)]
    pStop_(ngettext(length(missing),"parameter","parameters")," ",
      paste(sQuote(missing),collapse=","),
      " not found in ",sQuote("params"))
  }

  pompLoad(object)

  ker <- reuman.kernel(kernel.width)
  discrep <- spect.discrep(object,ker=ker,weights=weights)

  ofun <- function (par) {
    params[idx] <- par
    coef(object,transform=TRUE) <<- params
    object@simspec <- compute.spect.sim(object,vars=object@vars,
      params=object@params,nsim=object@nsim,seed=object@seed,
      transform.data=object@transform.data,detrend=object@detrend,ker=ker)
    discrep <<- spect.discrep(object,ker=ker,weights=weights)
    if (is.finite(discrep) || is.na(fail.value)) discrep else fail.value
  }

  environment(ofun) <- list2env(
    list(object=object,fail.value=fail.value,
      params=params,idx=idx,discrep=discrep,seed=seed,ker=ker,
      weights=weights),
    parent=parent.frame(2)
  )

  new("spect_match_objfun",ofun,env=environment(ofun))
}

## compute a measure of the discrepancies between simulations and data
spect.discrep <- function (object, ker, weights) {

  discrep <- array(dim=c(length(object@freq),length(object@vars)))
  sim.means <- colMeans(object@simspec)
  for (j in seq_along(object@freq)) {
    for (k in seq_along(object@vars)) {
      discrep[j,k] <- ((object@datspec[j,k]-sim.means[j,k])^2)/
        mean((object@simspec[,j,k]-sim.means[j,k])^2)
    }
    discrep[j,] <- weights[j]*discrep[j,]
  }

  sum(discrep)

}

##' @name spect-spect_match_objfun
##' @rdname spect
##' @aliases spect,spect_match_objfun-method
##' @export
setMethod(
  "spect",
  signature=signature(data="spect_match_objfun"),
  definition=function (data, ..., verbose=getOption("verbose", FALSE)) {
    spect(data@env$object,...,seed=data@env$seed,verbose=verbose)
  }
)

##' @name logLik-spect_match_objfun
##' @rdname loglik
##' @aliases logLik,spect_match_objfun-method
##' @export
setMethod(
  "logLik",
  signature=signature(object="spect_match_objfun"),
  definition=function (object) {
    discrep <- spect.discrep(object$env$object,ker=object@env$ker,
      weights=object@env$weights)
    -discrep
  }
)
