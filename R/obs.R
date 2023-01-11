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
##' @param format format of the returned object
##' @param \dots ignored
##' @export
setMethod(
  "obs",
  signature=signature(object="pomp"),
  definition=function (object, vars, ...,
    format = c("array", "data.frame")) {
    varnames <- rownames(object@data)
    if (missing(vars))
      vars <- varnames
    else if (!all(vars%in%varnames))
      pStop("obs","some elements of ",
        sQuote("vars")," correspond to no observed variable.")
    y <- object@data[vars,,drop=FALSE]
    dimnames(y) <- setNames(list(vars,NULL),c("name",object@timename))
    format <- match.arg(format)
    if (format == "data.frame") {
      y <- data.frame(time=time(object),t(y))
      names(y)[1L] <- object@timename
    }
    y
  }
)

##' @rdname obs
##' @importFrom dplyr bind_rows
##' @export
setMethod(
  "obs",
  signature=signature(object="listie"),
  definition=function (object, vars, ...,
    format = c("array", "data.frame")) {
    format <- match.arg(format)
    y <- lapply(object,obs,vars=vars,format=format,...)
    if (format == "data.frame") {
      bind_rows(y,.id=".id")
    } else {
      y
    }
  }
)
