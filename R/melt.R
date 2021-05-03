##' @rdname melt
##' @keywords internal
##' @docType import
##' @importFrom reshape2 melt
##' @details See \code{\link[reshape2:melt]{reshape2::melt}} for details.
##' @export
reshape2::melt

##' @rdname melt
##' @name melt,pomp-method
##' @keywords internal
##' @details A \sQuote{pomp} object can be melted into a data frame.
##' @inheritParams reshape2::melt
##' @export
setMethod(
  "melt",
  signature=signature(data="pomp"),
  definition=function (data, ...) {
    melt(as(data,"data.frame"),...)
  }
)
