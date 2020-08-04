##' Filtering mean
##'
##' The mean of the filtering distribution
##'
##' The filtering distribution is that of
##' \deqn{X(t_k) \vert Y(t_1)=y^*_1,\dots,Y(t_k)=y^*_k,}{Xk | Y1=y1*,\dots,Yk=yk*,}
##' where \eqn{X(t_k)}{Xk}, \eqn{Y(t_k)}{Yk} are the latent state and observable processes, respectively, and \eqn{y^*_t}{yt*} is the data, at time \eqn{t_k}{tk}.
##'
##' The filtering mean is therefore the expectation of this distribution
##' \deqn{E[X(t_k) \vert Y(t_1)=y^*_1,\dots,Y(t_k)=y^*_k].}{E[Xk | Y1=y1*,\dots,Yk=yk*].}
##'
##' @name filter.mean
##' @docType methods
##' @aliases filter.mean filter.mean,ANY-method filter.mean,missing-method
##' @include pfilter.R kalman.R
##' @rdname filter_mean
##' @family particle_filter_methods
##'
##' @param object result of a filtering computation
##' @param vars optional character; names of variables
##' @param \dots ignored
##'
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
    reqd_arg("filter.mean","object")
  }
)

setMethod(
  "filter.mean",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    undef_method("filter.mean",object)
  }
)

##' @name filter.mean-kalmand_pomp
##' @aliases filter.mean,kalmand_pomp-method
##' @rdname filter_mean
##' @export
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
##' @export
setMethod(
  "filter.mean",
  signature=signature(object="pfilterd_pomp"),
  definition=function (object, vars, ...) {
    if (missing(vars)) vars <- rownames(object@filter.mean)
    object@filter.mean[vars,,drop=FALSE]
  }
)
