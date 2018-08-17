##' Probe values
##'
##' Return the values of the probes, evaluated on simulations, and on the data.
##'
##' @name Probe values
##' @rdname probevals
##' @include probe.R
##' @aliases probevals probevals,missing-method probevals,ANY-method
##' @family summary statistics
NULL

setGeneric(
  "probevals",
  function (object, ...)
    standardGeneric("probevals")
)
setMethod(
  "probevals",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("probevals"),": ",sQuote("object"),
      " is a required argument",call.=FALSE)
  }
)

setMethod(
  "probevals",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("probevals")," is not defined for objects of class ",
      sQuote(class(object)),call.=FALSE)
  }
)

##' @name probevals-probed_pomp
##' @aliases probevals probevals,probed_pomp-method
##' @rdname probevals
##'
##' @param object an object of class \sQuote{probed_pomp}, or of a class that extends this.
##' @param \dots ignored
##'
setMethod(
  "probevals",
  signature=signature(object="probed_pomp"),
  definition=function (object, ...) {
    dv <- object@datvals
    list(
      data=array(dv,dim=c(1,length(dv)),
        dimnames=list(rep="data",probe=names(dv))),
      sims=object@simvals
    )
  }
)
