##' Forecast mean
##'
##' Mean of the one-step-ahead forecasting distribution.
##'
##' @name forecast
##' @rdname forecast
##' @aliases forecast forecast,missing-method forecast,ANY-method
##' @include kalman.R

setGeneric(
  "forecast",
  function (object, ...)
    standardGeneric("forecast")
)

setMethod(
  "forecast",
  signature=signature(object="missing"),
  definition=function (...) {
    reqd_arg("forecast","object")
  }
)

setMethod(
  "forecast",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    undef_method("forecast",object)
  }
)

##' @name forecast-kalmand_pomp
##' @aliases forecast,kalmand_pomp-method
##' @rdname forecast
##' @inheritParams filter.mean-kalmand_pomp
##'
##' @export
setMethod(
  "forecast",
  signature=signature(object="kalmand_pomp"),
  definition=function (object, vars, ...) {
    if (missing(vars)) vars <- rownames(object@forecast)
    object@forecast[vars,,drop=FALSE]
  }
)
