##' pomp_stop
##'
##' Custom error function
##' @name pomp_stop
##' @keywords internal
##'
##' @param fn name of function
##' @param \dots message
##'
pomp_stop <- function (fn = NULL, ...) {
  if (!is.null(fn))
    stop("in ",sQuote(fn),": ",...,call.=FALSE)
  else
    stop(...,call.=FALSE)
}
