##' Undefined
##'
##' Check for undefined methods.
##'
##' @return
##' Returns \code{TRUE} if the \pkg{pomp} workhorse method is undefined,
##' \code{FALSE} if it is defined,
##' and \code{NA} if the question is inapplicable.
##'
##' @name undefined
##' @rdname undefined
##' @include pstop.R
##'
##' @param object  object to test.
##'
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
