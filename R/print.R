##' Print methods
##'
##' These methods print their argument and return it *invisibly*.
##'
##' @name print
##' @rdname print
##' @include show.R
##' @aliases print,listie-method print,pomp_fun-method print,unshowable-method
NULL

setGeneric(
    "print",
    function (x, ...)
        standardGeneric("print")
)

##' @export
setMethod(
  "print",
  signature=signature(x="unshowable"),
  definition=function (x, ...) {
    show(x)
    invisible(x)
  }
)

##' @export
setMethod(
  "print",
  signature=signature(x="listie"),
  definition=function (x, ...) {
    show(x)
    invisible(x)
  }
)

##' @export
setMethod(
  "print",
  "pomp_fun",
  function (x, ...) show(x)
)
