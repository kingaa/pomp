##' Latent states
##'
##' Extract the latent states from a \sQuote{pomp} object.
##'
##' @name states
##' @docType methods
##' @include pomp_class.R
##' @rdname states
NULL

setGeneric(
  "states",
  function (object, ...)
    standardGeneric("states")
)

##' @name states-pomp
##' @aliases states states,pomp-method
##' @rdname states
##' @inheritParams obs
##'
setMethod(
  "states",
  signature=signature(object="pomp"),
  definition=function (object, vars, ...) {
    if (length(object@states)==0) {
      NULL
    } else {
      varnames <- rownames(object@states)
      if (missing(vars)) vars <- varnames
      else if (!all(vars%in%varnames))
        stop("in ",sQuote("states"),": some elements of ",
          sQuote("vars")," correspond to no state variable",call.=FALSE)
      x <- object@states[vars,,drop=FALSE]
      dimnames(x) <- setNames(list(vars,time(object)),
        c("variable",object@timename))
      x
    }
  }
)
