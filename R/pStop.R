##' pStop
##'
##' Custom error function
##' @name pStop
##' @keywords internal
##'
##' @param fn name of function
##' @param \dots message
##'
pStop <- function (fn, ...) {
  fn <- as.character(fn)
  if (length(fn) > 0 && nzchar(fn[1L]))
    stop("in ",sQuote(fn[1L]),": ",...,call.=FALSE)
  else
    stop(...,call.=FALSE)
}

pStop_ <- function (...) {
    stop(...,call.=FALSE)
}

pWarn <- function (fn, ...) {
  fn <- as.character(fn)
  if (length(fn) > 0 && nzchar(fn[1L]))
    warning("in ",sQuote(fn[1L]),": ",...,call.=FALSE)
  else
    warning(...,call.=FALSE)
}

pWarn_ <- function (...) {
  warning(...,call.=FALSE)
}

