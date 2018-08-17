##' Filtering mean
##'
##' The mean of the filtering distribution
##'
##' The filtering distribution is that of
##' \deqn{X_t \vert y_1,\dots,y_t,}{X_t | y_1,\dots,y_t,}
##' where \eqn{X_t}, \eqn{y_t} are the state vector and data, respectively,
##' at time \eqn{t}. 
##'
##' @name filter.mean
##' @aliases filter.mean filter.mean,ANY-method filter.mean,missing-method
##' @include pfilter.R kalman.R
##' @rdname filter_mean
##' @family particle filter methods
NULL

setGeneric(
  "filter.mean",
  function (object, ...)
    standardGeneric("filter.mean")
)

setMethod(
  "filter.mean",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("filter.mean"),": ",sQuote("object"),
      " is a required argument",call.=FALSE)
  }
)

setMethod(
  "filter.mean",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("filter.mean")," is not defined for objects of class ",
      sQuote(class(object)),call.=FALSE)
  }
)

##' @name filter.mean-kalmand_pomp
##' @aliases filter.mean,kalmand_pomp-method
##' @rdname filter_mean
##'
##' @param object an object of class \sQuote{pomp}, or of a class extending \sQuote{pomp}
##' @param vars optional character; names of variables
##' @param \dots ignored
##'
setMethod(
  "filter.mean",
  signature=signature(object="kalmand_pomp"),
  definition=function (object, vars, ...) {
    if (missing(vars)) vars <- rownames(object@filter.mean)
    object@filter.mean[vars,,drop=FALSE]
  }
)

##' @name filter.mean-pfilterd_pomp
##' @aliases filter.mean,pfilterd_pomp-method
##' @rdname filter_mean
setMethod(
  "filter.mean",
  signature=signature(object="pfilterd_pomp"),
  definition=function (object, vars, ...) {
    if (missing(vars)) vars <- rownames(object@filter.mean)
    object@filter.mean[vars,,drop=FALSE]
  }
)
