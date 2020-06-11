##' @rdname melt
##' @keywords internal
##' @importFrom reshape2 melt
##' 
##' @details See \code{reshape2::\link[reshape2]{melt}} for details.
##' 
##' @export
reshape2::melt

##' @rdname melt
##' @name melt,pomp-method
##' @keywords internal
##' @details A \sQuote{pomp} object can be melted into a data frame.
##' @inheritParams reshape2::melt
##' 
##' @export
setMethod(
  "melt",
  signature=signature(data="pomp"),
  definition=function (data, ...) {
    melt(as(data,"data.frame"),...)
  }
)
