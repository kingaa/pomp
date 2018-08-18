##' Conditional log likelihood
##'
##' The estimated conditional log likelihood from a fitted model.
##'
##' The conditional log likelihood is defined to be the value of the density
##' of \deqn{Y_t | Y_1,\dots,Y_{t-1}}{Yt | Y1,\dots,Y(t-1)} evaluated at \eqn{Y_t = y^*_t}{Yt = yt*}.
##' Here, \eqn{Y_t}{Yt} is the observable process and \eqn{y^*_t}{yt*} is the data, at time \eqn{t}.
##'
##' \deqn{\ell_t(\theta) = \mathrm{Prob}[y_t \vert y_1, \dots, y_{t-1}],}{ell_t(theta)=Prob[y_t | y_1, \dots, y_(t-1)],} where \eqn{y_t} are the data, at time \eqn{t}.
##'
##' @name cond.logLik
##' @docType methods
##' @rdname cond_logLik
##' @include pomp_class.R kalman.R pfilter.R
##' @aliases cond.logLik cond.logLik,missing-method cond.logLik,ANY-method
##' @family particle filter methods
##'
##' @return
##' numerical value of the conditional log likelihood.
##' Note that some methods compute not the log likelihood itself but instead a related quantity.
##' To keep the code simple, the \code{logLik} function is nevertheless used to extract this quantity.
NULL

setGeneric("cond.logLik",
  function (object, ...)
    standardGeneric("cond.logLik")
)

setMethod(
  "cond.logLik",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("cond.logLik"),": ",sQuote("object"),
      " is a required argument",call.=FALSE)
  }
)

setMethod(
  "cond.logLik",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("cond.logLik")," is not defined for objects of class ",
      sQuote(class(object)),call.=FALSE)
  }
)

##' @name cond.logLik-kalmand_pomp
##' @aliases cond.logLik,kalmand_pomp-method
##' @rdname cond_logLik
##'
##' @param object an object of class \sQuote{pomp}, or of a class extending \sQuote{pomp}
##' @param \dots ignored
##'
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
