## extract prediction variances

setGeneric(
    "pred.var",
    function (object, ...)
        standardGeneric("pred.var")
)

setMethod(
  "pred.var",
  signature=signature(object="pfilterd_pomp"),
  definition=function (object, vars, ...) {
    if (missing(vars)) vars <- rownames(object@pred.var)
    object@pred.var[vars,,drop=FALSE]
  }
)

setMethod(
  "pred.var",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("pred.var")," is not defined for objects of class ",
      sQuote(class(object)),call.=FALSE)
  }
)
