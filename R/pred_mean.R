##' Prediction mean
##'
##' The mean of the prediction distribution
##'
##' The prediction distribution is that of
##' \deqn{X(t_k) \vert Y(t_1)=y^*_1,\dots,Y(t_{k-1})=y^*_{k-1},}{Xk | Y1=y1*,\dots,Y(k-1)=y(k-1)*,}
##' where \eqn{X(t_k)}{Xk}, \eqn{Y(t_k)}{Yk} are the latent state and observable processes, respectively, and \eqn{y^*_k}{yk*} is the data, at time \eqn{t_k}{tk}.
##'
##' The prediction mean is therefore the expectation of this distribution
##' \deqn{E[X(t_k) \vert Y(t_1)=y^*_1,\dots,Y(t_{k-1})=y^*_{k-1}].}{E[Xk | Y1=y1*,\dots,Y(k-1)=y(k-1)*].}
##'
##' @name pred.mean
##' @aliases pred.mean pred.mean,ANY-method pred.mean,missing-method
##' @include pfilter.R kalman.R
##' @rdname pred_mean
##' @family particle_filter_methods
##' @inheritParams filter.mean
##'
NULL

setGeneric(
  "pred.mean",
  function (object, ...)
    standardGeneric("pred.mean")
)

setMethod(
  "pred.mean",
  signature=signature(object="missing"),
  definition=function (...) {
    reqd_arg("pred.mean","object")
  }
)

setMethod(
  "pred.mean",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    undef_method("pred.mean",object)
  }
)

##' @name pred.mean-kalmand_pomp
##' @aliases pred.mean,kalmand_pomp-method
##' @rdname pred_mean
##' @export
setMethod(
  "pred.mean",
  signature=signature(object="kalmand_pomp"),
  definition=function (object, vars, ...) {
    if (missing(vars)) vars <- rownames(object@pred.mean)
    object@pred.mean[vars,,drop=FALSE]
  }
)

##' @name pred.mean-pfilterd_pomp
##' @aliases pred.mean,pfilterd_pomp-method
##' @rdname pred_mean
##' @export
setMethod(
  "pred.mean",
  signature=signature(object="pfilterd_pomp"),
  definition=function (object, vars, ...) {
    if (missing(vars)) vars <- rownames(object@pred.mean)
    object@pred.mean[vars,,drop=FALSE]
  }
)
