##' Log likelihood
##'
##' Extract the estimated log likelihood (or related quantity) from a fitted model.
##'
##' @name logLik
##' @rdname loglik
##' @aliases logLik logLik,ANY-method logLik,missing-method
##' logLik,listie-method
##' @include pfilter.R mif2.R pmcmc.R probe.R kalman.R nlf.R listie.R
##' @include objfun.R spect_match.R nlf.R
##'
##' @param object fitted model object
##' @param \dots ignored
##'
##' @return numerical value of the log likelihood.
##' Note that some methods compute not the log likelihood itself but instead a related quantity.
##' To keep the code simple, the \code{logLik} function is nevertheless used to extract this quantity.
##'
NULL

##' @rdname loglik
##' @export
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

##' @rdname loglik
##' @export
setMethod(
  "logLik",
  signature=signature(object="listie"),
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

##' @name logLik-probed_pomp
##' @aliases logLik,probed_pomp-method
##' @rdname loglik
##'
##' @return
##' When \code{object} is of \sQuote{probed_pomp} class (i.e., the result of a \code{probe} computation), \code{logLik} retrieves the \dQuote{synthetic likelihood} (see \code{\link{probe}}).
##'
setMethod(
  "logLik",
  signature=signature(object="probed_pomp"),
  definition=function(object)object@synth.loglik
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

##' @name logLik-objfun
##' @rdname loglik
##' @aliases logLik,objfun-method
##' @export
setMethod(
  "logLik",
  signature=signature(object="objfun"),
  definition=function (object) {
    object@env$loglik
  }
)

##' @name logLik-spect_match_objfun
##' @rdname loglik
##' @aliases logLik,spect_match_objfun-method
##' @export
setMethod(
  "logLik",
  signature=signature(object="spect_match_objfun"),
  definition=function (object) {
    -object@env$discrep
  }
)

##' @name logLik-nlf_objfun
##' @aliases logLik,nlf_objfun-method
##' @rdname loglik
##'
##' @return
##' When \code{object} is an NLF objective function, i.e., the result of a call to \code{nlf_objfun},
##' \code{logLik} retrieves the \dQuote{quasi log likelihood} (see \code{\link{nlf}}).
##'
##' @export
setMethod(
  "logLik",
  signature=signature(object="nlf_objfun"),
  definition = function(object, ...) {
    object@env$logql
  }
)
