##' obs
##'
##' Extract the data array from a \sQuote{pomp} object.
##'
##' @name obs
##' @aliases obs,ANY-method obs,missing-method
##' @docType methods
##' @rdname obs
##' @include pomp_class.R
##' @importFrom stats setNames
##' @family extraction methods
##'
NULL

setGeneric(
    "obs",
    function (object, ...)
        standardGeneric("obs")
)

setMethod(
  "obs",
  signature=signature(object="missing"),
  definition=function (object, ...) {
    reqd_arg("obs","object")
  }
)

setMethod(
  "obs",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    undef_method("obs","object")
  }
)

##' @rdname obs
##' @param object an object of class \sQuote{pomp}, or of a class extending \sQuote{pomp}
##' @param vars names of variables to retrieve
##' @param \dots ignored
##' @export
setMethod(
  "obs",
  signature=signature(object="pomp"),
  definition=function (object, vars, ...) {
    varnames <- rownames(object@data)
    if (missing(vars))
      vars <- varnames
    else if (!all(vars%in%varnames))
      pStop("obs","some elements of ",
        sQuote("vars")," correspond to no observed variable.")
    y <- object@data[vars,,drop=FALSE]
    dimnames(y) <- setNames(list(vars,NULL),c("variable","time"))
    y
  }
)

##' @rdname obs
##' @export
setMethod(
  "obs",
  signature=signature(object="listie"),
  definition=function (object, vars, ...) {
    lapply(object,obs,vars=vars)
  }
)
