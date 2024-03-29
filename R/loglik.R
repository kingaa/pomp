##' Log likelihood
##'
##' Extract the estimated log likelihood (or related quantity) from a fitted model.
##'
##' @name logLik
##' @rdname loglik
##' @aliases logLik,ANY-method logLik,missing-method
##' @include pfilter.R wpfilter.R mif2.R pmcmc.R probe.R kalman.R nlf.R listie.R
##' @include objfun.R spect_match.R nlf.R
##' @family extraction methods
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
##' @importFrom stats logLik
##' @export
setGeneric("logLik")

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

##' @rdname loglik
##' @return
##' When \code{object} is of \sQuote{pfilterd_pomp} class (i.e., the result of a \code{wpfilter} computation), \code{logLik} retrieves the estimated log likelihood.
##' @export
setMethod(
  "logLik",
  signature=signature(object="pfilterd_pomp"),
  definition=function(object)object@loglik
)

##' @rdname loglik
##' @return
##' When \code{object} is of \sQuote{wpfilterd_pomp} class (i.e., the result of a \code{wpfilter} computation), \code{logLik} retrieves the estimated log likelihood.
##' @export
setMethod(
  "logLik",
  signature=signature(object="wpfilterd_pomp"),
  definition=function(object)object@loglik
)

##' @rdname loglik
##' @return
##' When \code{object} is of \sQuote{probed_pomp} class (i.e., the result of a \code{\link{probe}} computation), \code{logLik} retrieves the \dQuote{synthetic likelihood}.
##' @export
setMethod(
  "logLik",
  signature=signature(object="probed_pomp"),
  definition=function(object)object@synth.loglik
)

##' @rdname loglik
##' @return
##' When \code{object} is of \sQuote{kalmand_pomp} class (i.e., the result of an \code{\link{eakf}} or \code{\link{enkf}} computation), \code{logLik} retrieves the estimated log likelihood.
##' @export
setMethod(
  "logLik",
  signature=signature(object="kalmand_pomp"),
  definition=function(object)object@loglik
)

##' @rdname loglik
##' @return
##' When \code{object} is of \sQuote{pmcmcd_pomp} class (i.e., the result of a \code{\link{pmcmc}} computation), \code{logLik} retrieves the estimated log likelihood as of the last particle filter operation.
##' @export
setMethod(
  "logLik",
  signature=signature(object="pmcmcd_pomp"),
  definition=function (object)
    object@loglik
)

##' @rdname loglik
##' @return
##' When \code{object} is of \sQuote{bsmcd_pomp} class (i.e., the result of a \code{\link{bsmc2}} computation), \code{logLik} retrieves the \dQuote{log evidence}.
##' @export
setMethod(
  "logLik",
  signature=signature(object="bsmcd_pomp"),
  definition=function(object)object@log.evidence
)

##' @rdname loglik
##' @export
setMethod(
  "logLik",
  signature=signature(object="objfun"),
  definition=function (object) {
    object@env$loglik
  }
)

##' @rdname loglik
##' @return
##' When \code{object} is of \sQuote{spect_match_objfun} class (i.e., an objective function constructed by \code{\link{spect_objfun}}), \code{logLik} retrieves minus the spectrum mismatch.
##' @export
setMethod(
  "logLik",
  signature=signature(object="spect_match_objfun"),
  definition=function (object) {
    -object@env$discrep
  }
)

##' @rdname loglik
##' @return
##' When \code{object} is an NLF objective function, i.e., the result of a call to \code{\link{nlf_objfun}},
##' \code{logLik} retrieves the \dQuote{quasi log likelihood}.
##' @export
setMethod(
  "logLik",
  signature=signature(object="nlf_objfun"),
  definition = function(object, ...) {
    object@env$logql
  }
)
