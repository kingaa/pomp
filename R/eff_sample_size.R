##' Effective sample size
##'
##' Estimate the effective sample size of a Monte Carlo computation.
##'
##' Effective sample size is computed as
##' \deqn{\left(\sum_i\!w_{it}^2\right)^{-1},}{1/(sum(w_it^2)),}
##' where \eqn{w_{it}}{w_it} is the normalized weight of particle \eqn{i} at time \eqn{t}.
##'
##' @name eff_sample_size
##' @rdname eff_sample_size
##' @include pfilter.R wpfilter.R bsmc2.R
##' @aliases eff_sample_size,missing-method eff_sample_size,ANY-method
##' @family particle filter methods
##' @family extraction methods
##' @inheritParams filter_mean
##'
NULL

setGeneric(
  "eff_sample_size",
  function (object, ...)
    standardGeneric("eff_sample_size")
)

setMethod(
  "eff_sample_size",
  signature=signature(object="missing"),
  definition=function (...) {
    reqd_arg("eff_sample_size","object")
  }
)

setMethod(
  "eff_sample_size",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    undef_method("eff_sample_size",object)
  }
)

##' @rdname eff_sample_size
##' @export
setMethod(
  "eff_sample_size",
  signature=signature(object="bsmcd_pomp"),
  definition=function (object, ...,
    format = c("numeric", "data.frame")) {
    format <- match.arg(format)
    if (format == "numeric") {
      object@eff.sample.size
    } else {
      data.frame(
        time=time(object),
        eff.sample.size=object@eff.sample.size
      )
    }
  }
)

##' @rdname eff_sample_size
##' @export
setMethod(
  "eff_sample_size",
  signature=signature(object="pfilterd_pomp"),
  definition=function (object, ...,
    format = c("numeric", "data.frame")) {
    format <- match.arg(format)
    if (format == "numeric") {
      object@eff.sample.size
    } else {
      data.frame(
        time=time(object),
        eff.sample.size=object@eff.sample.size
      )
    }
  }
)

##' @rdname eff_sample_size
##' @export
setMethod(
  "eff_sample_size",
  signature=signature(object="wpfilterd_pomp"),
  definition=function (object, ...,
    format = c("numeric", "data.frame")) {
    format <- match.arg(format)
    if (format == "numeric") {
      object@eff.sample.size
    } else {
      data.frame(
        time=time(object),
        eff.sample.size=object@eff.sample.size
      )
    }
  }
)

##' @rdname eff_sample_size
##' @export
setMethod(
  "eff_sample_size",
  signature=signature(object="pfilterList"),
  definition=function (object, ...,
    format = c("numeric", "data.frame")) {
    format <- match.arg(format)
    x <- lapply(object,eff_sample_size,format=format)
    if (format == "data.frame") {
      x <- rbind_fill(x,.id=".id")
    }
    x
  }
)
