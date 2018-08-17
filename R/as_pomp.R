##' as.pomp
##'
##' Coerce to a \sQuote{pomp} object
##'
##' @name as.pomp
##' @rdname as_pomp
##' @include pomp_class.R
##' @keywords internal
##'
##' @param object the object to be coerced
##' @param \dots additional arguments

as.pomp <- function (object, ...) {
  as(object,"pomp")
}
