##' Prediction mean
##'
##' The mean of the prediction distribution
##'
##' The prediction distribution is that of
##' \deqn{X_t \vert y_1,\dots,y_{t-1},}{X_t | y_1,\dots,y_(t-1),}
##' where \eqn{X_t}, \eqn{y_t} are the state vector and data, respectively,
##' at time \eqn{t}.
##' 
##' @name pred.mean
##' @aliases pred.mean pred.mean,ANY-method pred.mean,missing-method
##' @include pfilter.R kalman.R
##' @rdname pred_mean
##' @family particle filter methods
NULL

setGeneric(
  "pred.mean",
  function (object, ...)
    standardGeneric("pred.mean")
)

setMethod(
  "pred.mean",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("pred.mean"),": ",sQuote("object"),
      " is a required argument",call.=FALSE)
  }
)

setMethod(
  "pred.mean",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("pred.mean")," is not defined for objects of class ",
      sQuote(class(object)),call.=FALSE)
  }
)

##' @name pred.mean-kalmand_pomp
##' @aliases pred.mean,kalmand_pomp-method
##' @rdname pred_mean
##' @inheritParams filter.mean-kalmand_pomp
##'
setMethod(
  "pred.mean",
  signature=signature(object="kalmand_pomp"),
  definition=function (object, vars, ...) {
    if (missing(vars)) vars <- rownames(object@pred.mean)
    object@pred.mean[vars,,drop=FALSE]
  }
)

##' @name pred.mean-pfilterd_pomp
##' @aliases pred.mean,pfilterd_pomp-method
##' @rdname pred_mean
##' @inheritParams filter.mean-kalmand_pomp
##'
setMethod(
  "pred.mean",
  signature=signature(object="pfilterd_pomp"),
  definition=function (object, vars, ...) {
    if (missing(vars)) vars <- rownames(object@pred.mean)
    object@pred.mean[vars,,drop=FALSE]
  }
)
