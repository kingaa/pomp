## a simple method to extract the data array
setGeneric(
    "obs",
    function (object, ...)
        standardGeneric("obs")
)

setMethod(
  "obs",
  signature=signature(object="pomp"),
  definition=function (object, vars, ...) {
    varnames <- rownames(object@data)
    if (missing(vars))
      vars <- varnames
    else if (!all(vars%in%varnames))
      stop("in ",sQuote("obs"),": some elements of ",
        sQuote("vars")," correspond to no observed variable",call.=FALSE)
    y <- object@data[vars,,drop=FALSE]
    dimnames(y) <- setNames(list(vars,time(object)),
      c("variable",object@timename))
    y
  }
)
