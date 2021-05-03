##' The zero time
##'
##' Get and set the zero-time.
##'
##' @name timezero
##' @rdname timezero
##' @docType methods
##' @aliases timezero<- timezero,missing-method timezero,ANY-method
##' timezero<-,missing-method timezero<-,ANY-method
##' @return
##' the value of the zero time
##'
NULL

setGeneric(
  "timezero",
  function (object, ...)
    standardGeneric("timezero")
)

setGeneric(
  "timezero<-",
  function (object, ..., value)
    standardGeneric("timezero<-")
)

##' @rdname timezero
##' @param object an object of class \sQuote{pomp}, or of a class that extends \sQuote{pomp}
##' @param \dots ignored
##' @export
setMethod(
  "timezero",
  signature=signature(object="pomp"),
  definition = function (object, ...) object@t0
)

##' @rdname timezero
##' @param value numeric; the new zero-time value
##' @export
setMethod(
  "timezero<-",
  signature=signature(object="pomp"),
  definition=function(object,...,value) {
    ep <- "timezero<-"
    if (!(is.numeric(value) && length(value) == 1L && is.finite(value)))
      pStop(ep,"the zero-time ",sQuote("t0")," must be a single finite number.")
    if (value > object@times[1L])
      pStop(ep,"the zero-time ",sQuote("t0"),
        " must occur no later than the first observation.")
    storage.mode(value) <- "double"
    object@t0 <- value
    object
  }
)
