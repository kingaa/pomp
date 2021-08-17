##' Print methods
##'
##' These methods print their argument and return it *invisibly*.
##'
##' @name print
##' @rdname print
##' @keywords internal
##' @param x object to print
##' @param ... ignored
##' @include show.R
NULL

setGeneric("print")

##' @rdname print
##' @export
setMethod(
  "print",
  signature=signature(x="unshowable"),
  definition=function (x, ...) {
    show(x)
    invisible(x)
  }
)

##' @rdname print
##' @export
setMethod(
  "print",
  signature=signature(x="listie"),
  definition=function (x, ...) {
    show(x)
    invisible(x)
  }
)

##' @rdname print
##' @export
setMethod(
  "print",
  "pomp_fun",
  function (x, ...) {
    show(x)
    invisible(x)
  }
)
