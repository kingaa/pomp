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
  if (!is.null(fn) || !nzchar(fn[1L]))
    stop("in ",sQuote(fn[1L]),": ",...,call.=FALSE)
  else
    stop(...,call.=FALSE)
}

pomp_warn <- function (fn = NULL, ...) {
  if (!is.null(fn) || !nzchar(fn[1L]))
    warning("in ",sQuote(fn[1L]),": ",...,call.=FALSE)
  else
    warning(...,call.=FALSE)
}
