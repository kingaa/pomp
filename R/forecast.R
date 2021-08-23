##' Forecast mean
##'
##' Mean of the one-step-ahead forecasting distribution.
##'
##' @name forecast
##' @rdname forecast
##' @aliases forecast,missing-method forecast,ANY-method
##' @family extraction methods
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

##' @rdname forecast
##' @inheritParams filter.mean
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
