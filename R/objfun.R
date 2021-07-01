##' Objective functions
##'
##' Methods common to \pkg{pomp} stateful objective functions
##'
##' @name objfun
##' @rdname objfun
##' @keywords internal
##' @include traj_match.R spect_match.R probe_match.R nlf.R
##' @include loglik.R summary.R coef.R as_data_frame.R as_pomp.R
##'
##' @section Important Note:
##' Since \pkg{pomp} cannot guarantee that the \emph{final} call an optimizer makes to the function is a call \emph{at} the optimum, it cannot guarantee that the parameters stored in the function are the optimal ones.
##' Therefore, it is a good idea to evaluate the function on the parameters returned by the optimization routine, which will ensure that these parameters are stored.
##'
NULL

setClassUnion(
  "objfun",
  members=c(
    "traj_match_objfun",
    "spect_match_objfun",
    "probe_match_objfun",
    "nlf_objfun"
  )
)

##' @rdname coef
##' @export
setMethod(
  "coef",
  signature=signature(object="objfun"),
  definition=function (object, ...) {
    coef(object@env$object,...)
  }
)

##' @rdname summary
##' @export
setMethod(
  "summary",
  signature=signature(object="objfun"),
  definition=function (object, ...) {
    summary(object@env$object)
  }
)

setAs(
  from="objfun",
  to="pomp",
  def = function (from) {
    as(from@env$object,"pomp")
  }
)

setAs(
  from="objfun",
  to="data.frame",
  def = function (from) {
    as(as(from@env$object,"pomp"),"data.frame")
  }
)

##' @rdname partrans
##' @export
setMethod(
  "partrans",
  signature=signature(object="objfun"),
  definition=function (object, ...) {
    partrans(as(object,"pomp"),...)
  }
)

##' @rdname simulate
##' @export
setMethod(
  "simulate",
  signature=signature(object="objfun"),
  definition=function (object, seed = NULL, ...) {
    simulate(as(object,"pomp"),seed=seed,...)
  }
)

##' @rdname probe
##' @export
setMethod(
  "probe",
  signature=signature(data="objfun"),
  definition=function (data, seed = NULL, ...) {
    probe(as(data,"pomp"),seed=seed,...)
  }
)

##' @rdname pfilter
##' @export
setMethod(
  "pfilter",
  signature=signature(data="objfun"),
  definition=function (data, ...) {
    pfilter(as(data,"pomp"),...)
  }
)

##' @rdname spect
##' @export
setMethod(
  "spect",
  signature=signature(data="objfun"),
  definition=function (data, seed = NULL, ...) {
    spect(as(data,"pomp"),seed=seed,...)
  }
)
