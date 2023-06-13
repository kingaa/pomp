##' Latent states
##'
##' Extract the latent states from a \sQuote{pomp} object.
##'
##' @name states
##' @aliases states,ANY-method states,missing-method
##' @rdname states
##' @docType methods
##' @include pomp_class.R melt.R
##' @importFrom stats setNames
##' @family extraction methods
##'
NULL

setGeneric(
  "states",
  function (object, ...)
    standardGeneric("states")
)

setMethod(
  "states",
  signature=signature(object="missing"),
  definition=function (object, ...) {
    reqd_arg("states","object")
  }
)

setMethod(
  "states",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    undef_method("states","object")
  }
)

##' @rdname states
##' @inheritParams obs
##' @export
setMethod(
  "states",
  signature=signature(object="pomp"),
  definition=function (object, vars, ...,
    format = c("array", "data.frame")) {
    if (length(object@states)==0) {
      x <- array(NA_real_,dim=c(0,length(object@times)),
        dimnames=setNames(list(NULL,NULL),c("name",object@timename)))
    } else if (missing(vars)) {
      x <- object@states
      names(dimnames(x)) <- c("name",object@timename)
    } else {
      varnames <- rownames(object@states)
      if (!all(vars%in%varnames))
        pStop(who="states","some elements of ",
          sQuote("vars")," correspond to no state variable")
      x <- object@states[vars,,drop=FALSE]
      dimnames(x) <- setNames(list(vars,NULL),c("name",object@timename))
    }
    format <- match.arg(format)
    if (format == "data.frame") {
      x <- data.frame(time=time(object),t(x))
      names(x)[1L] <- object@timename
    }
    x
  }
)

##' @rdname states
##' @inheritParams obs
##' @export
setMethod(
  "states",
  signature=signature(object="listie"),
  definition=function (object, vars, ...,
    format = c("array", "data.frame")) {
    format <- match.arg(format)
    x <- lapply(object,states,vars=vars,format=format,...)
    if (format == "data.frame") {
      rbind_fill(x,.id=".id")
    } else {
      x
    }
  }
)
