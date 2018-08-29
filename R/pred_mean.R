##' Prediction mean
##'
##' The mean of the prediction distribution
##'
##' The prediction distribution is that of
##' \deqn{X_t \vert Y_1=y^*_1,\dots,Y_{t-1}=y^*_{t-1},}{Xt | Y1=y1*,\dots,Y(t-1)=y(t-1)*,}
##' where \eqn{X_t}{Xt}, \eqn{Y_t}{Yt} are the latent state and observable processes, respectively, and \eqn{y^*_t}{yt*} is the data, at time \eqn{t}.
##'
##' The prediction mean is therefore the expectation of this distribution
##' \deqn{E[X_t \vert Y_1=y^*_1,\dots,Y_{t-1}=y^*_{t-1}].}{E[Xt | Y1=y1*,\dots,Y(t-1)=y(t-1)*].}
##'
##' @name pred.mean
##' @aliases pred.mean pred.mean,ANY-method pred.mean,missing-method
##' @include pfilter.R kalman.R
##' @rdname pred_mean
##' @family particle filter methods
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
