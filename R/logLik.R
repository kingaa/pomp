## extract the estimated log likelihood

setGeneric(
    "logLik",
    function (object, ...)
        standardGeneric("logLik")
)

setMethod(
  "logLik",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("logLik"),": ",sQuote("object"),
      " is a required argument",call.=FALSE)
  }
)

setMethod(
  "logLik",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
      NA_real_
  }
)

setMethod(
  "logLik",
  signature=signature(object="listies"),
  definition=function(object, ...) {
    do.call(c,lapply(object,logLik))
  }
)

setMethod(
  "logLik",
  signature=signature(object="pfilterd_pomp"),
  definition=function(object,...)object@loglik
)

setMethod(
  "logLik",
  signature=signature(object="kalmand_pomp"),
  definition=function(object,...)object@loglik
)

setMethod(
  "logLik",
  signature=signature(object="pmcmcd_pomp"),
  definition=function (object, ...)
    object@loglik
)

setMethod(
  "logLik",
  signature=signature(object="probed_pomp"),
  definition=function(object,...)object@synth.loglik
)
