##' @importFrom reshape2 melt
##' @docType import
##' @export
reshape2::melt

##' @name melt,pomp-method
##' @include melt.R
##' @rdname as_data_frame
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
