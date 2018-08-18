##' Objective functions
##'
##' Methods common to \pkg{pomp} stateful objective functions
##'
##' @name objfun
##' @rdname objfun
##' @keywords internal
##' @include traj_match.R spect_match.R probe_match.R
##' @include loglik.R summary.R coef.R
##'
##' @section Important Note:
##' Since \pkg{pomp} cannot guarantee that the \emph{final} call an optimizer makes to the function is a call \emph{at} the optimum, it cannot guarantee that the parameters stored in the function are the optimal ones.
##' The best practice is to evaluate the function on the parameters returned by the optimization routine, which will ensure that these parameters are stored.
##'
NULL

setClassUnion(
  "objfun",
  members=c(
    "traj_match_objfun",
    "spect_match_objfun",
    "probe_match_objfun"
  )
)

##' @name coef-objfun
##' @rdname coef
##' @aliases coef,objfun-method
##' @export
setMethod(
  "coef",
  signature=signature(object="objfun"),
  definition=function (object, ...) {
    coef(object@env$object,...)
  }
)

##' @name summary-objfun
##' @rdname summary
##' @aliases summary,objfun-method
##' @export
setMethod(
  "summary",
  signature=signature(object="objfun"),
  definition=function (object) {
    summary(object@env$object)
  }
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

##' @name coerce-objfun-pomp
##' @aliases coerce,objfun,pomp-method
##' @rdname as_pomp
##'
setAs(
  from="objfun",
  to="pomp",
  def = function (from) {
    from@env$object
  }
)
