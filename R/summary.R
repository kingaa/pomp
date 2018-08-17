##' Summary methods
##'
##' Display a summary of a fitted model object.
##'
##' @rdname summary
##' @name summary
##'
##' @param object a fitted model object
##' @param \dots ignored
NULL

setGeneric(
  "summary",
  function (object, ...)
    standardGeneric("summary")
)
