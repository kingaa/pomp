##' Conjoin
##'
##' Combine two or more \sQuote{pomp} objects into a list-like \sQuote{listie}.
##'
##' @name conjoin
##' @rdname conjoin
##' @aliases c
##' @include listie.R concat.R
##' @param ... elements to be recursively combined into a \sQuote{listie}
##' @keywords internal
##' @details
##' \code{conjoin} applied to one or lists applies \code{c} to convert the list into a \sQuote{listie}.
##' In particular, \code{conjoin(A,B,C)} is equivalent to \code{do.call(c,unlist(list(A,B,C)))}.
##' 
NULL

##' @rdname conjoin
##' @export
c.Pomp <- function (...) concat(...)

##' @rdname conjoin
##' @export
conjoin <- function (...) {
  do.call(concat,unlist(list(...)))
}
