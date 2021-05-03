##' as.pomp
##'
##' Coerce to a \sQuote{pomp} object
##'
##' @name as.pomp
##' @docType methods
##' @keywords internal
##' @rdname as_pomp
##' @include pomp_class.R
##'
##' @param object the object to be coerced
##' @param \dots additional arguments
##' 
##' @export
as.pomp <- function (object, ...) {
  as(object,"pomp")
}
