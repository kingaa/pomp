##' Prediction variance
##'
##' The variance of the prediction distribution
##'
##' @name Prediction variance
##' @aliases pred.var pred.var,ANY-method pred.var,missing-method
##' @include pfilter.R kalman.R
##' @rdname pred_var
NULL

setGeneric(
  "pred.var",
  function (object, ...)
    standardGeneric("pred.var")
)

setMethod(
  "pred.var",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("pred.var"),": ",sQuote("object"),
      " is a required argument",call.=FALSE)
  }
)

setMethod(
  "pred.var",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("pred.var")," is not defined for objects of class ",
      sQuote(class(object)),call.=FALSE)
  }
)

##' @name pred.var-pfilterd_pomp
##' @aliases pred.var,pfilterd_pomp-method
##' @rdname pred_var
##' @inheritParams filter.mean-kalmand_pomp
##'
setMethod(
  "pred.var",
  signature=signature(object="pfilterd_pomp"),
  definition=function (object, vars, ...) {
    if (missing(vars)) vars <- rownames(object@pred.var)
    object@pred.var[vars,,drop=FALSE]
  }
)
