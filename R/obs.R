##' obs
##'
##' Extract the data array from a \sQuote{pomp} object.
##'
##' @name obs
##' @docType methods
##' @rdname obs
##' @include pomp_class.R
NULL

setGeneric(
    "obs",
    function (object, ...)
        standardGeneric("obs")
)

##' @name obs-pomp
##' @aliases obs obs,pomp-method
##' @rdname obs
##' @param object an object of class \sQuote{pomp}, or of a class extending \sQuote{pomp}
##' @param vars names of variables to retrieve
##' @param \dots ignored
##'
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
