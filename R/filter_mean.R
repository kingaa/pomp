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
##' @name filter_mean
##' @docType methods
##' @aliases filter_mean,ANY-method filter_mean,missing-method
##' @include pfilter.R kalman.R melt.R
##' @rdname filter_mean
##' @family particle filter methods
##' @family extraction methods
##'
##' @param object result of a filtering computation
##' @param vars optional character; names of variables
##' @param format format of the returned object.
##' @param ... ignored
##'
NULL

setGeneric(
  "filter_mean",
  function (object, ...)
    standardGeneric("filter_mean")
)

setMethod(
  "filter_mean",
  signature=signature(object="missing"),
  definition=function (...) {
    reqd_arg("filter_mean","object")
  }
)

setMethod(
  "filter_mean",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    undef_method("filter_mean",object)
  }
)

##' @rdname filter_mean
##' @export
setMethod(
  "filter_mean",
  signature=signature(object="kalmand_pomp"),
  definition=function (object, vars, ...,
    format = c("array", "data.frame")) {
    if (missing(vars)) vars <- rownames(object@filter.mean)
    format <- match.arg(format)
    if (format == "array") {
      object@filter.mean[vars,,drop=FALSE]
    } else {
      x <- melt(object@filter.mean[vars,,drop=FALSE])
      x$time <- time(object)[as.integer(x$time)]
      x
    }
  }
)

##' @rdname filter_mean
##' @export
setMethod(
  "filter_mean",
  signature=signature(object="pfilterd_pomp"),
  definition=function (object, vars, ...,
    format = c("array", "data.frame")) {
    if (missing(vars)) vars <- rownames(object@filter.mean)
    format <- match.arg(format)
    if (format == "array") {
      object@filter.mean[vars,,drop=FALSE]
    } else {
      x <- melt(object@filter.mean[vars,,drop=FALSE])
      x$time <- time(object)[as.integer(x$time)]
      x
    }
  }
)
