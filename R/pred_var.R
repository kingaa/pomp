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
##' @name pred.var
##' @aliases pred.var pred.var,ANY-method pred.var,missing-method
##' @include pfilter.R kalman.R
##' @rdname pred_var
##' @family particle filter methods
##' @inheritParams filter.mean-kalmand_pomp
##'
NULL

setGeneric(
  "pred.var",
  function (object, ...)
    standardGeneric("pred.var")
)

setMethod(
  "pred.var",
  signature=signature(object="missing"),
  definition=function (...) {
    reqd_arg("pred.var","object")
  }
)

setMethod(
  "pred.var",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    undef_method("pred.var",object)
  }
)

##' @name pred.var-pfilterd_pomp
##' @aliases pred.var,pfilterd_pomp-method
##' @rdname pred_var
##' @export
setMethod(
  "pred.var",
  signature=signature(object="pfilterd_pomp"),
  definition=function (object, vars, ...) {
    if (missing(vars)) vars <- rownames(object@pred.var)
    object@pred.var[vars,,drop=FALSE]
  }
)
