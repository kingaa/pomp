##' Deprecated functions
##'
##' @name deprecated
##' @rdname deprecated
##' @keywords internal
##' @include traces.R rprocess_spec.R as_data_frame.R
NULL

##' @name conv.rec
##' @rdname deprecated
##' @keywords internal
##' @export
conv.rec <- function (object, ...) {
  pWarn_(sQuote("conv.rec")," is deprecated and will be removed in a ",
    "forthcoming release.  Please use ",sQuote("traces")," instead.")
  traces(object,...)
}

##' @name values
##' @rdname deprecated
##' @keywords internal
##' @export
values <- function (object, ...) {
  pWarn_(sQuote("values")," is deprecated and will be removed in a ",
    "forthcoming release.  Please use ",sQuote("as.data.frame")," instead.")
  as(object,"data.frame")
}

##' @name onestep.dens
##' @rdname deprecated
##' @keywords internal
##' @export
onestep.dens <- function (dens.fun, PACKAGE) {
  pWarn_(sQuote("onestep.dens")," is deprecated and will be removed in a ",
    "forthcoming release.  ","Specify ",sQuote("dprocess")," directly ",
    "using a C snippet or R function.")
  if (!missing(PACKAGE))
    warning("in ",sQuote("onestep.dens"),": ",sQuote("PACKAGE")," ignored.",
      call.=FALSE)
  dens.fun
}

##' @name onestep.sim
##' @rdname deprecated
##' @keywords internal
##' @export
onestep.sim <- function (step.fun) {
  pWarn_(sQuote("onestep.sim")," is deprecated and will be removed in a ",
    "forthcoming release.  Use ",sQuote("onestep")," instead.")
  onestep(step.fun)
}

##' @name discrete.time.sim
##' @rdname deprecated
##' @keywords internal
##' @export
discrete.time.sim <- function (step.fun, delta.t = 1) {
  pWarn_(sQuote("discrete.time.sim")," is deprecated and will be removed in a ",
    "forthcoming release.  Use ",sQuote("discrete_time")," instead.")
  discrete_time(step.fun=step.fun,delta.t=delta.t)
}

##' @name euler.sim
##' @rdname deprecated
##' @keywords internal
##' @export
euler.sim <- function (step.fun, delta.t) {
  pWarn_(sQuote("euler.sim")," is deprecated and will be removed in a ",
    "forthcoming release.  Use ",sQuote("euler")," instead.")
  euler(step.fun=step.fun,delta.t=delta.t)
}

##' @name gillespie.sim
##' @rdname deprecated
##' @keywords internal
##' @export
gillespie.sim <- function (rate.fun, v, hmax = Inf) {
  pWarn_(sQuote("gillespie.sim")," is deprecated and will be removed in a ",
    "forthcoming release.  Use ",sQuote("gillespie")," instead.")
  gillespie(rate.fun=rate.fun,v=v,hmax=hmax)
}

##' @name gillespie.hl.sim
##' @rdname deprecated
##' @keywords internal
##' @export
gillespie.hl.sim <- function (..., .pre = "", .post = "", hmax = Inf) {
  pWarn_(sQuote("gillespie.hl.sim")," is deprecated and will be removed in a ",
    "forthcoming release.  Use ",sQuote("gillespie_hl")," instead.")
  gillespie_hl(...,.pre=.pre,.post=.post,hmax=hmax)
}
