##' Sampling the prior and posterior distributions
##'
##' These functions return the samples of the prior and posterior distributions
##' from a \code{bsmc2} run.
##'
##' @name samples
##' @rdname samples
##' @include pomp-package.R bsmc2.R
##' @aliases prior_samples prior_samples,missing-method prior_samples,ANY-method
##' posterior_samples posterior_samples,missing-method posterior_samples,ANY-method
NULL

setGeneric("prior_samples",
  function(object,...)standardGeneric("prior_samples"))

setGeneric("posterior_samples",
  function(object,...)standardGeneric("posterior_samples"))

setMethod(
  "prior_samples",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("prior_samples"),": ",sQuote("object"),
      " is a required argument",call.=FALSE)
  }
)

setMethod(
  "prior_samples",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("prior_samples")," is not defined for objects of class ",
      sQuote(class(object)),call.=FALSE)
  }
)

setMethod(
  "posterior_samples",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("posterior_samples"),": ",sQuote("object"),
      " is a required argument",call.=FALSE)
  }
)

setMethod(
  "posterior_samples",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("posterior_samples")," is not defined for objects of class ",
      sQuote(class(object)),call.=FALSE)
  }
)

##' @name prior_samples-bsmcd_pomp
##' @aliases prior_samples prior_samples,bsmcd_pomp-method
##' @rdname samples
##' @param object the result of a \code{bsmc2} computation
##' @param pars names of parameters
##' @param \dots ignored
##'
setMethod(
  "prior_samples",
  signature=signature(object="bsmcd_pomp"),
  definition=function (object, pars, ...) {
    if (missing(pars)) pars <- object@est
    object@prior[pars,,drop=FALSE]
  }
)

##' @name posterior_samples-bsmcd_pomp
##' @aliases posterior_samples,bsmcd_pomp-method
##' @rdname samples
setMethod(
  "posterior_samples",
  signature=signature(object="bsmcd_pomp"),
  definition=function (object, pars, ...) {
    if (missing(pars)) pars <- object@est
    object@post[pars,,drop=FALSE]
  }
)
