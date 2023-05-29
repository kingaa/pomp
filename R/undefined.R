##' Undefined
##'
##' Check for undefined methods.
##'
##' @name undefined
##' @rdname undefined
##' @return
##' Returns \code{TRUE} if the \pkg{pomp} workhorse method is undefined,
##' \code{FALSE} if it is defined,
##' and \code{NA} if the question is inapplicable.
##' @include pstop.R
##' @keywords internal
##' @param object  object to test.
##' @param ... currently ignored.
##'
NULL

setGeneric(
  "undefined",
  function (object, ...)
    standardGeneric("undefined")
)

setMethod(
  "undefined",
  signature=signature(object="NULL"),
  definition=function (object, ...) TRUE
)

setMethod(
  "undefined",
  signature=signature(object="ANY"),
  definition=function (object, ...) NA
)

setMethod(
  "undefined",
  signature=signature(object="missing"),
  definition=function (...) NA
)
