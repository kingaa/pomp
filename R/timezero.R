##' The zero time
##'
##' Get and set the zero-time.
##'
##' @name timezero
##' @rdname timezero
##' @docType methods
##' @aliases timezero timezero<- timezero,missing-method timezero,ANY-method
##' timezero<-,missing-method timezero<-,ANY-method
##'
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

##' @name timezero-pomp
##' @aliases timezero,pomp-method
##' @rdname timezero
##'
##' @param object an object of class \sQuote{pomp}, or of a class that extends \sQuote{pomp}
##' @param \dots ignored
##'
setMethod(
  "timezero",
  signature=signature(object="pomp"),
  definition=function(object,...)object@t0
)

##' @name timezero<--pomp
##' @aliases timezero<-,pomp-method
##' @rdname timezero
##'
##' @param value numeric; the new zero-time value
##'
setMethod(
  "timezero<-",
  signature=signature(object="pomp"),
  definition=function(object,...,value) {
    ep <- paste0("in ",sQuote("timezero<-"),": ")
    if (!(is.numeric(value) && length(value) == 1L && is.finite(value)))
      stop(ep,"the zero-time ",sQuote("t0"),
        " must be a single finite number.",call.=FALSE)
    if (value > object@times[1L])
      stop(ep,"the zero-time ",sQuote("t0"),
        " must occur no later than the first observation.",call.=FALSE)
    storage.mode(value) <- "double"
    object@t0 <- value
    object
  }
)
