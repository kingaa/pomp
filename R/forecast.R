##' Forecast mean
##'
##' Mean of the one-step-ahead forecasting distribution.
##'
##' @name Forecast
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
    stop("in ",sQuote("forecast"),": ",sQuote("object"),
      " is a required argument",call.=FALSE)
  }
)

setMethod(
  "forecast",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("forecast")," is not defined for objects of class ",
      sQuote(class(object)),call.=FALSE)
  }
)

##' @name forecast-kalmand_pomp
##' @aliases forecast,kalmand_pomp-method
##' @rdname forecast
##' @inheritParams filter.mean-kalmand_pomp
##'

setMethod(
  "forecast",
  signature=signature(object="kalmand_pomp"),
  definition=function (object, vars, ...) {
    if (missing(vars)) vars <- rownames(object@forecast)
    object@forecast[vars,,drop=FALSE]
  }
)
