setGeneric(
    "eff.sample.size",
    function (object, ...)
        standardGeneric("eff.sample.size")
)

setMethod(
  "eff.sample.size",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("eff.sample.size"),": ",sQuote("object"),
      " is a required argument",call.=FALSE)
  }
)

setMethod(
  "eff.sample.size",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("eff.sample.size")," is not defined for objects of class ",
      sQuote(class(object)),call.=FALSE)
  }
)

setMethod(
  "eff.sample.size",
  signature=signature(object="bsmcd_pomp"),
  definition=function(object,...)object@eff.sample.size
)

setMethod(
  "eff.sample.size",
  signature=signature(object="pfilterd_pomp"),
  definition=function(object,...)object@eff.sample.size
)
