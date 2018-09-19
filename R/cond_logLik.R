##' Conditional log likelihood
##'
##' The estimated conditional log likelihood from a fitted model.
##'
##' The conditional likelihood is defined to be the value of the density
##' of \deqn{Y_t | Y_1,\dots,Y_{t-1}}{Yt | Y1,\dots,Y(t-1)} evaluated at \eqn{Y_t = y^*_t}{Yt = yt*}.
##' Here, \eqn{Y_t}{Yt} is the observable process and \eqn{y^*_t}{yt*} is the data, at time \eqn{t}.
##'
##' Thus the conditional log likelihood at time \eqn{t} is
##' \deqn{\ell_t(\theta) = \log f[Y_t=y^*_t \vert Y_1=y^*_1, \dots, Y_{t-1}=y^*_{t-1}],}{ell_t(theta)=log f[Yt = yt*t | Y1=y1*, \dots, Y(t-1)=y(t-1)*],}
##' where \eqn{f} is the probability density above.
##'
##' @name cond.logLik
##' @docType methods
##' @rdname cond_logLik
##' @include pomp_class.R kalman.R pfilter.R
##' @aliases cond.logLik cond.logLik,missing-method cond.logLik,ANY-method
##' @family particle filter methods
##' @inheritParams filter.mean
##'
##' @return
##' The numerical value of the conditional log likelihood.
##' Note that some methods compute not the log likelihood itself but instead a related quantity.
##' To keep the code simple, the \code{cond.logLik} function is nevertheless used to extract this quantity.
NULL

setGeneric("cond.logLik",
  function (object, ...)
    standardGeneric("cond.logLik")
)

setMethod(
  "cond.logLik",
  signature=signature(object="missing"),
  definition=function (...) {
    reqd_arg("cond.logLik","object")
  }
)

setMethod(
  "cond.logLik",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    undef_method("cond.logLik",object)
  }
)

##' @name cond.logLik-kalmand_pomp
##' @aliases cond.logLik,kalmand_pomp-method
##' @rdname cond_logLik
##' @export
setMethod(
  "cond.logLik",
  signature=signature(object="kalmand_pomp"),
  definition=function(object,...)object@cond.loglik
)

##' @name cond.logLik-pfilterd_pomp
##' @aliases cond.logLik,pfilterd_pomp-method
##' @rdname cond_logLik
##' @export
setMethod(
  "cond.logLik",
  signature=signature(object="pfilterd_pomp"),
  definition=function(object,...)object@cond.loglik
)

##' @name cond.logLik-bsmcd_pomp
##' @aliases cond.logLik,bsmcd_pomp-method
##' @rdname cond_logLik
##'
##' @return
##' When \code{object} is of class \sQuote{bsmcd_pomp} (i.e., the result of a \code{bsmc2} computation), \code{cond.logLik} returns the conditional log \dQuote{evidence} (see \code{\link{bsmc2}}).
##'
##' @export
setMethod(
  "cond.logLik",
  signature=signature(object="bsmcd_pomp"),
  definition=function(object,...)object@cond.log.evidence
)
