setGeneric(
    "filter.mean",
    function (object, ...)
        standardGeneric("filter.mean")
)

setMethod(
  "filter.mean",
  signature=signature(object="kalmand_pomp"),
  definition=function (object, vars, ...) {
    if (missing(vars)) vars <- rownames(object@filter.mean)
    object@filter.mean[vars,,drop=FALSE]
  }
)

setMethod(
  "filter.mean",
  signature=signature(object="pfilterd_pomp"),
  definition=function (object, vars, ...) {
    if (missing(vars)) vars <- rownames(object@filter.mean)
    object@filter.mean[vars,,drop=FALSE]
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
