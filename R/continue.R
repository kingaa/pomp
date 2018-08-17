##' Continue an iterative calculation
##'
##' Continue an iterative computation where it left off.
##'
##' @name continue
##' @aliases continue continue,missing-method continue,ANY-method
##' @rdname continue
##'
##' @param object the result of an iterative \pkg{pomp} computation
##' @param \dots additional arguments will be passed to the underlying method.
##' This allows one to modify parameters used in the original computations.
NULL

##' @name continue
##' @rdname continue
setGeneric(
  "continue",
  function (object, ...)
    standardGeneric("continue")
)

setMethod(
  "continue",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("continue"),": ",sQuote("object"),
      " is a required argument",call.=FALSE)
  }
)

setMethod(
  "continue",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("continue")," is not defined for objects of class ",
      sQuote(class(object)),call.=FALSE)
  }
)
