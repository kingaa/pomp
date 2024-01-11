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
##' @section Warning! Objective functions based on C snippets:
##' If you use C snippets (see \code{\link{Csnippet}}), a dynamically loadable library will be built.
##' As a rule, \pkg{pomp} functions load this library as needed and unload it when it is no longer needed.
##' The stateful objective functions are an exception to this rule.
##' For efficiency, calls to the objective function do not execute \code{\link{pompLoad}} or \code{\link{pompUnload}}:
##' rather, it is assumed that \code{\link{pompLoad}} has been called before any call to the objective function.
##' When a stateful objective function using one or more C snippets is created, \code{\link{pompLoad}} is called internally to build and load the library:
##' therefore, within a single \R session, if one creates a stateful objective function, one can freely call that objective function and (more to the point) pass it to an optimizer that calls it freely, without needing to call \code{\link{pompLoad}}.
##' On the other hand, if one retrieves a stored objective function from a file, or passes one to another \R session, one must call \code{\link{pompLoad}} before using it.
##' \strong{Failure to do this will typically result in a segmentation fault (i.e., it will crash the \R session).}
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

##' @rdname coef
##' @export
setMethod(
  "coef<-",
  signature=signature(object="objfun"),
  definition=function (object, pars, transform = FALSE, ..., value) {
    coef(object@env$object,pars=pars,transform=transform,...) <- value
    object@env$params <- coef(object@env$object,transform=TRUE)
    object
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

setIs(
  class1="objfun",
  class2="pomp",
  coerce = function (from) {
    as(from@env$object,"pomp")
  },
  replace= function (from, to, value) {
    pStop(who=-2L,"cannot replace the pomp object in a stateful objective function.")
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

##' @rdname load
##' @export
setMethod(
  "pompLoad",
  signature=signature(object="objfun"),
  definition = function (object, ...) {
    object@env$.gnsi <- TRUE
    pompLoad_internal(object@env$object,...)
  })

##' @rdname load
##' @export
setMethod(
  "pompUnload",
  signature=signature(object="objfun"),
  definition = function (object, ...) {
    object@env$.gnsi <- TRUE
    pompUnload_internal(object@env$object,...)
  })
