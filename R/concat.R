##' Concatenate
##'
##' Concatenate two or more \sQuote{pomp} objects into a list-like \sQuote{listie}.
##'
##' @name concat
##' @rdname concat
##' @aliases c
##' @include listie.R conc.R
##' @param ... elements to be recursively combined into a \sQuote{listie}
##' @details
##' \code{concat} applied to one or more \sQuote{pomp} objects or lists of \sQuote{pomp} objects converts the list into a \sQuote{listie}.
##' In particular, \code{concat(A,B,C)} is equivalent to \code{do.call(c,unlist(list(A,B,C)))}.
##' 
##' @example examples/concat.R
##' 
NULL

##' @rdname concat
##' @export
c.Pomp <- function (...) conc(...)

##' @rdname concat
##' @export
concat <- function (...) {
  do.call(conc,unlist(list(...)))
}
