##' Combine
##'
##' Combine two or more \sQuote{pomp} objects into a list-like \sQuote{listie}.
##'
##' @name combine
##' @rdname combine
##' @aliases c
##' @include listie.R concat.R
##' @param ... elements to be recursively combined into a \sQuote{listie}
##' @keywords internal
##' @return
##' \code{combine} applied to one or lists applies \code{c} to convert the list into a \sQuote{listie}.
##' In particular, \code{combine(A,B,C)} is equivalent to \code{do.call(c,list(A,B,C))}.
##' 
NULL

##' @rdname combine
##' @export
c.Pomp <- function (...) concat(...)

##' @rdname combine
##' @export
combine <- function (...) {
  do.call(concat,list(...))
}
