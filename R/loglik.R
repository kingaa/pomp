##' Log likelihood
##'
##' Extract the estimated log likelihood from a fitted model.
##'
##' @name logLik
##' @rdname loglik
##' @aliases logLik logLik,ANY-method logLik,missing-method
##' logLik,listies-method
##' @include pfilter.R mif2.R pmcmc.R probe.R kalman.R nlf.R listies.R
##'
##' @param object fitted model object
##' @param \dots ignored
##'
##' @return numerical value of the log likelihood.
##' Note that some methods compute not the log likelihood itself but instead a related quantity.
##' To keep the code simple, the \code{logLik} function is nevertheless used to extract this quantity.
##'
NULL

##' @name logLik-generic
##' @rdname loglik
setGeneric(
  "logLik",
  function (object, ...)
    standardGeneric("logLik")
)

setMethod(
  "logLik",
  signature=signature(object="missing"),
  definition=function (...) {
    reqd_arg("logLik","object")
  }
)

setMethod(
  "logLik",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    NA_real_
  }
)

##' @export
setMethod(
  "logLik",
  signature=signature(object="listies"),
  definition=function(object,...) {
    do.call(c,lapply(object,logLik,...))
  }
)

##' @name logLik-pfilterd_pomp
##' @aliases logLik,pfilterd_pomp-method
##' @rdname loglik
##' @export
setMethod(
  "logLik",
  signature=signature(object="pfilterd_pomp"),
  definition=function(object)object@loglik
)

##' @name logLik-kalmand_pomp
##' @aliases logLik,kalmand_pomp-method
##' @rdname loglik
##' @export
setMethod(
  "logLik",
  signature=signature(object="kalmand_pomp"),
  definition=function(object)object@loglik
)

##' @name logLik-pmcmcd_pomp
##' @aliases logLik,pmcmcd_pomp-method
##' @rdname loglik
##' @export
setMethod(
  "logLik",
  signature=signature(object="pmcmcd_pomp"),
  definition=function (object)
    object@loglik
)

##' @name logLik-bsmcd_pomp
##' @aliases logLik,bsmcd_pomp-method
##' @rdname loglik
##'
##' @importFrom stats logLik
##'
##' @return
##' When \code{object} is of \sQuote{bsmcd_pomp} class (i.e., the result of a \code{bsmc2} computation), \code{logLik} retrieves the \dQuote{log evidence} (see \code{\link{bsmc2}}).
##'
##' @export
setMethod(
  "logLik",
  signature=signature(object="bsmcd_pomp"),
  definition=function(object)object@log.evidence
)
