##' Melt
##'
##' Convert arrays, lists, and other objects to data frames.
##'
##' \code{melt} converts its first argument to a data frame.
##' It is a simplified version of the \code{melt} command provided by the now-deprecated \pkg{reshape2} package.
##'
##' @name melt
##' @rdname melt
##' @keywords internal
##' @aliases melt,ANY-method melt,missing-method
##' @param data object to convert
##' @param ... ignored
##'
NULL

setGeneric(
  "melt",
  function (data, ...)
    standardGeneric("melt")
)

setMethod(
  "melt",
  signature=signature(data="missing"),
  definition=function (data, ...) {
    reqd_arg("melt","data")
  }
)

##' @rdname melt
##' @importFrom stats setNames
##' @export
setMethod(
  "melt",
  signature=signature(data="ANY"),
  definition=function (data, ...) {
    data.frame(value=data)
  }
)

##' @rdname melt
##' @details An array can be melted into a data frame.
##' In this case, the data frame will have one row per entry in the array.
##' @export
setMethod(
  "melt",
  signature=signature(data="array"),
  definition=function (data, ...) {
    dn <- dimnames(data)
    nullnames <- which(unlist(lapply(dn,is.null)))
    dn[nullnames] <- lapply(nullnames,function(i)seq_len(dim(data)[i]))
    labels <- expand.grid(dn,KEEP.OUT.ATTRS=FALSE,stringsAsFactors=FALSE)
    dim(data) <- NULL
    cbind(labels,value=data)
  }
)

##' @rdname melt
##' @details A list can be melted into a data frame.
##' This operation is recursive.
##' A variable will be appended to distinguish the separate list entries.
##' @importFrom dplyr bind_rows
##' @export
setMethod(
  "melt",
  signature=signature(data="list"),
  definition=function (data, ..., level = 1) {
    if (length(unique(lapply(data,mode))) > 1L)
      pStop("melt","refusing to melt data of incompatible types.")
    bind_rows(
      lapply(data,melt,level=level+1,...),
      .id=paste0("L",level)
    )
  }
)
