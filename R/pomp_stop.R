##' pomp_stop
##'
##' Custom error function
##' @name pomp_stop
##' @keywords internal
NULL

pomp_stop <- function (..., which = 1) {
  f <- sQuote(as.character(sys.call(which)[[1]]))
  stop("in ",f,": ",...,call.=FALSE)
}
