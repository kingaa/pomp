##' Effective sample size
##'
##' Estimate the effective sample size of a Monte Carlo computation.
##'
##' Effective sample size is computed as
##' \deqn{\left(\sum_i\!w_{it}^2\right)^{-1},}{1/(sum(w_it^2)),}
##' where \eqn{w_{it}}{w_it} is the normalized weight of particle \eqn{i} at time \eqn{t}.
##'
##' @name eff.sample.size
##' @rdname eff_sample_size
##' @include pfilter.R wpfilter.R bsmc2.R
##' @aliases eff.sample.size,missing-method eff.sample.size,ANY-method
##' @family particle filter methods
##' @inheritParams filter.mean
##'
NULL

setGeneric(
    "eff.sample.size",
    function (object, ...)
        standardGeneric("eff.sample.size")
)

setMethod(
  "eff.sample.size",
  signature=signature(object="missing"),
  definition=function (...) {
    reqd_arg("eff.sample.size","object")
  }
)

setMethod(
  "eff.sample.size",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    undef_method("eff.sample.size",object)
  }
)

##' @rdname eff_sample_size
##' @export
setMethod(
  "eff.sample.size",
  signature=signature(object="bsmcd_pomp"),
  definition=function(object,...)object@eff.sample.size
)

##' @rdname eff_sample_size
##' @export
setMethod(
  "eff.sample.size",
  signature=signature(object="pfilterd_pomp"),
  definition=function(object,...)object@eff.sample.size
)

##' @rdname eff_sample_size
##' @export
setMethod(
  "eff.sample.size",
  signature=signature(object="wpfilterd_pomp"),
  definition=function(object,...)object@eff.sample.size
)
