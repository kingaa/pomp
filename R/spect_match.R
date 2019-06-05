##' Spectrum matching
##'
##' Estimation of parameters by matching power spectra
##'
##' In spectrum matching, one attempts to minimize the discrepancy between a \acronym{POMP} model's predictions and data, as measured in the frequency domain by the power spectrum.
##'
##' \code{spect_objfun} constructs an objective function that measures the discrepancy.
##' It can be passed to any one of a variety of numerical optimization routines, which will adjust model parameters to minimize the discrepancies between the power spectrum of model simulations and that of the data.
##'
##' @name spect.match
##' @docType methods
##' @rdname spect_match
##' @family pomp parameter estimation methods
##' @aliases spect_objfun spect_objfun,missing-method spect_objfun,ANY-method
##' @example examples/spect_match.R
##'
##' @include spect.R probe_match.R loglik.R plot.R
##'
##' @param weights optional numeric or function.
##' The mismatch between model and data is measured by a weighted average of mismatch at each frequency.
##' By default, all frequencies are weighted equally.
##' \code{weights} can be specified either as a vector (which must have length equal to the number of frequencies) or as a function of frequency.
##' If the latter, \code{weights(freq)} must return a nonnegative weight for each frequency.
##'
##' @inheritParams probe.match
##' @inheritParams spect
##' @inheritParams pomp
##'
##' @return
##' \code{spect_objfun} constructs a stateful objective function for spectrum matching.
##' Specifically, \code{spect_objfun} returns an object of class \sQuote{spect_match_objfun}, which is a function suitable for use in an \code{\link{optim}}-like optimizer.
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
    env="environment",
    est="character"
  )
)

setGeneric(
  "spect_objfun",
  function (data, ...)
    standardGeneric("spect_objfun")
)

setMethod(
  "spect_objfun",
  signature=signature(data="missing"),
  definition=function (...) {
    reqd_arg("spect_objfun","data")
  }
)

setMethod(
  "spect_objfun",
  signature=signature(data="ANY"),
  definition=function (data, ...) {
    undef_method("spect_objfun",data)
  }
)

##' @name spect_objfun-data.frame
##' @aliases spect_objfun,data.frame-method
##' @rdname spect_match
##' @export
setMethod(
  "spect_objfun",
  signature=signature(data="data.frame"),
  definition=function(data,
    est = character(0), weights = 1, fail.value = NA,
    vars, kernel.width, nsim, seed = NULL, transform.data = identity,
    detrend = c("none","mean","linear","quadratic"),
    params, rinit, rprocess, rmeasure, partrans,
    ..., verbose = getOption("verbose", FALSE)) {

    tryCatch(
      smof.internal(
        data,
        est=est,
        weights=weights,
        fail.value=fail.value,
        vars=vars,
        kernel.width=kernel.width,
        nsim=nsim,
        seed=seed,
        transform.data=transform.data,
        detrend=detrend,
        params=params,
        rinit=rinit,
        rprocess=rprocess,
        rmeasure=rmeasure,
        partrans=partrans,
        ...,
        verbose=verbose
      ),
      error = function (e) pStop("spect_objfun",conditionMessage(e))
    )

  }
)

##' @name spect_objfun-pomp
##' @aliases spect_objfun,pomp-method
##' @rdname spect_match
##' @export
setMethod(
  "spect_objfun",
  signature=signature(data="pomp"),
  definition=function(data,
    est = character(0), weights = 1, fail.value = NA,
    vars, kernel.width, nsim, seed = NULL, transform.data = identity,
    detrend = c("none","mean","linear","quadratic"),
    ..., verbose = getOption("verbose", FALSE)) {

    tryCatch(
      smof.internal(
        data,
        est=est,
        weights=weights,
        fail.value=fail.value,
        vars=vars,
        kernel.width=kernel.width,
        nsim=nsim,
        seed=seed,
        transform.data=transform.data,
        detrend=detrend,
        ...,
        verbose=verbose
      ),
      error = function (e) pStop("spect_objfun",conditionMessage(e))
    )

  }
)

##' @name spect_objfun-spectd_pomp
##' @aliases spect_objfun,spectd_pomp-method
##' @rdname spect_match
##' @export
setMethod(
  "spect_objfun",
  signature=signature(data="spectd_pomp"),
  definition=function(data,
    est = character(0), weights = 1, fail.value = NA,
    vars, kernel.width, nsim, seed = NULL, transform.data = identity,
    detrend,
    ..., verbose = getOption("verbose", FALSE)) {

    if (missing(vars)) vars <- data@vars
    if (missing(kernel.width)) kernel.width <- data@kernel.width
    if (missing(nsim)) nsim <- data@nsim
    if (missing(transform.data)) transform.data <- data@transform.data
    if (missing(detrend)) detrend <- data@detrend

    spect_objfun(
      as(data,"pomp"),
      est=est,
      weights=weights,
      fail.value=fail.value,
      vars=vars,
      kernel.width=kernel.width,
      nsim=nsim,
      seed=seed,
      transform.data=transform.data,
      detrend=detrend,
      ...,
      verbose=verbose
    )

  }
)

##' @name spect_objfun-spect_match_objfun
##' @aliases spect_objfun,spect_match_objfun-method
##' @rdname spect_match
##' @export
setMethod(
  "spect_objfun",
  signature=signature(data="spect_match_objfun"),
  definition=function(data,
    est, weights, fail.value, seed = NULL,
    ..., verbose = getOption("verbose", FALSE)) {

    if (missing(est)) est <- data@est
    if (missing(weights)) weights <-data@env$weights
    if (missing(fail.value)) fail.value <- data@env$fail.value

    spect_objfun(
      data@env$object,
      est=est,
      weights=weights,
      fail.value=fail.value,
      seed=seed,
      ...,
      verbose=verbose
    )

  }
)

smof.internal <- function (object,
  est, weights, fail.value,
  vars, kernel.width, nsim, seed, transform.data, detrend,
  ..., verbose) {

  verbose <- as.logical(verbose)

  object <- spect(object,vars=vars,kernel.width=kernel.width,
    nsim=nsim,seed=seed, transform.data=transform.data,detrend=detrend,
    ...,verbose=verbose)

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
      error = function (e)
        pStop_(sQuote("weights")," function: ",conditionMessage(e))
    )
  } else {
    pStop_(sQuote("weights"),
      " must be specified as a vector or as a function")
  }

  if (any(!is.finite(weights) | weights<0))
    pStop_(sQuote("weights")," should be nonnegative and finite")
  weights <- weights/mean(weights)

  params <- coef(object,transform=TRUE)

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
    object@simspec <- compute.spect.sim(
      object,
      vars=object@vars,
      params=object@params,
      nsim=object@nsim,
      seed=object@seed,
      transform.data=object@transform.data,
      detrend=object@detrend,
      ker=ker
    )
    discrep <<- spect.discrep(object,ker=ker,weights=weights)
    if (is.finite(discrep) || is.na(fail.value)) discrep else fail.value
  }

  environment(ofun) <- list2env(
    list(object=object,fail.value=fail.value,
      params=params,idx=idx,discrep=discrep,seed=seed,ker=ker,
      weights=weights),
    parent=parent.frame(2)
  )

  new("spect_match_objfun",ofun,env=environment(ofun),est=est)

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
##'
##' @details
##' When \code{spect} operates on a spectrum-matching objective function (a \sQuote{spect_match_objfun} object), by default, the
##' random-number generator seed is fixed at the value given when the objective function was constructed.
##' Specifying \code{NULL} or an integer for \code{seed} overrides this behavior.
##'
##' @export
setMethod(
  "spect",
  signature=signature(data="spect_match_objfun"),
  definition=function (data, seed,
    ..., verbose=getOption("verbose", FALSE)) {

    if (missing(seed)) seed <- data@env$seed

    spect(
      data@env$object,
      seed=seed,
      ...,
      verbose=verbose
    )

  }
)

##' @name coerce-spect_match_objfun-spectd_pomp
##' @aliases coerce,spect_match_objfun,spectd_pomp-method
##' @rdname spect
##'
setAs(
  from="spect_match_objfun",
  to="spectd_pomp",
  def = function (from) {
    from@env$object
  }
)

##' @name plot-spect_match_objfun
##' @aliases plot,spect_match_objfun-method
##' @rdname plot
##' @export
setMethod(
  "plot",
  signature=signature(x="spect_match_objfun"),
  definition=function (x, ...) {
    plot(as(x,"spectd_pomp"),...)
  }
)
