##' Prediction variance
##'
##' The variance of the prediction distribution
##'
##' The prediction distribution is that of
##' \deqn{X(t_k) \vert Y(t_1)=y^*_1,\dots,Y(t_{k-1})=y^*_{k-1},}{Xk | Y1=y1*,\dots,Y(k-1)=y(k-1)*,}
##' where \eqn{X(t_k)}{Xk}, \eqn{Y(t_k)}{Yk} are the latent state and observable processes, respectively, and \eqn{y^*_k}{yk*} is the data, at time \eqn{t_k}{tk}.
##'
##' The prediction variance is therefore the variance of this distribution
##' \deqn{\mathrm{Var}[X(t_k) \vert Y(t_1)=y^*_1,\dots,Y(t_{k-1})=y^*_{k-1}].}{Var[Xk | Y1=y1*,\dots,Y(k-1)=y(k-1)*].}
##'
##' @name pred_var
##' @aliases pred_var,ANY-method pred_var,missing-method
##' @include pfilter.R kalman.R
##' @rdname pred_var
##' @family particle filter methods
##' @family extraction methods
##' @inheritParams filter_mean
##'
NULL

setGeneric(
  "pred_var",
  function (object, ...)
    standardGeneric("pred_var")
)

setMethod(
  "pred_var",
  signature=signature(object="missing"),
  definition=function (...) {
    reqd_arg("pred_var","object")
  }
)

setMethod(
  "pred_var",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    undef_method("pred_var",object)
  }
)

##' @rdname pred_var
##' @export
setMethod(
  "pred_var",
  signature=signature(object="pfilterd_pomp"),
  definition=function (object, vars, ...,
    format = c("array", "data.frame")) {
    if (missing(vars)) vars <- rownames(object@pred.var)
    format <- match.arg(format)
    if (format == "array") {
      object@pred.var[vars,,drop=FALSE]
    } else {
      x <- melt(object@pred.var[vars,,drop=FALSE])
      x$time <- time(object)[as.integer(x$time)]
      x
    }
  }
)
