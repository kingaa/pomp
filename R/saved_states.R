##' Saved states
##'
##' Retrieve latent state trajectories from a particle filter calculation.
##'
##' When one calls \code{\link{pfilter}} with \code{save.states=TRUE}, the latent state vector associated with each particle is saved.
##' This can be extracted by calling \code{saved.states} on the \sQuote{pfilterd.pomp} object.
##' 
##' @name saved.states
##' @aliases saved.states saved.states,ANY-method saved.states,missing-method
##' @include pfilter.R pmcmc.R
##' @rdname saved_states
##' @family particle filter methods
##' @inheritParams filter.mean
##'
##' @return The saved states are returned in the form of a list, with one element per time-point.
##' Each element consists of a matrix, with one row for each state variable and one column for each particle.
##' 
NULL

setGeneric(
  "saved.states",
  function (object,...) standardGeneric("saved.states")
)

setMethod(
  "saved.states",
  signature=signature(object="missing"),
  definition=function (...) {
    reqd_arg("saved.states","object")
  }
)

setMethod(
  "saved.states",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    undef_method("saved.states",object)
  }
)

##' @name saved.states-pfilterd_pomp
##' @aliases saved.states,pfilterd_pomp-method
##' @rdname saved_states
##'
##' @export
setMethod(
  "saved.states",
  signature=signature(object="pfilterd_pomp"),
  definition=function (object, ...) {
    object@saved.states
  }
)

##' @name saved.states-pfilterList
##' @aliases saved.states,pfilterList-method
##' @rdname saved_states
##' @export
setMethod(
  "saved.states",
  signature=signature(object="pfilterList"),
  definition=function (object, ...) {
    lapply(object,saved.states,...)
  }
)
