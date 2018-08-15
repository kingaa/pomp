## extract prediction means

setGeneric(
    "pred.mean",
    function (object, ...)
        standardGeneric("pred.mean")
)

setMethod(
  "pred.mean",
  signature=signature(object="kalmand_pomp"),
  definition=function (object, vars, ...) {
    if (missing(vars)) vars <- rownames(object@pred.mean)
    object@pred.mean[vars,,drop=FALSE]
  }
)

setMethod(
  "pred.mean",
  signature=signature(object="pfilterd_pomp"),
  definition=function (object, vars, ...) {
    if (missing(vars)) vars <- rownames(object@pred.mean)
    object@pred.mean[vars,,drop=FALSE]
  }
)

setMethod(
  "pred.mean",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("pred.mean")," is not defined for objects of class ",
      sQuote(class(object)),call.=FALSE)
  }
)

