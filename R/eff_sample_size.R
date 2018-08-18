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
##' @include pfilter.R bsmc2.R
##' @aliases eff.sample.size eff.sample.size,missing-method eff.sample.size,ANY-method
##' @family particle filter methods
##'
##' @param object object extending class \sQuote{pomp}
##' @param \dots ignored
##'
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

##' @name eff.sample.size-bsmcd_pomp
##' @aliases eff.sample.size,bsmcd_pomp-method
##' @rdname eff_sample_size
##' @export
setMethod(
  "eff.sample.size",
  signature=signature(object="bsmcd_pomp"),
  definition=function(object,...)object@eff.sample.size
)

##' @name eff.sample.size-pfilterd_pomp
##' @aliases eff.sample.size,pfilterd_pomp-method
##' @rdname eff_sample_size
##' @export
setMethod(
  "eff.sample.size",
  signature=signature(object="pfilterd_pomp"),
  definition=function(object,...)object@eff.sample.size
)
