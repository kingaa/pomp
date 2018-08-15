setGeneric(
    "print",
    function (x, ...)
        standardGeneric("print")
)

setMethod(
  "print",
  signature=signature(x="unshowable"),
  definition=function (x, ...) {
    show(x)
    invisible(x)
  }
)

setMethod(
  "print",
  signature=signature(x="listies"),
  definition=function (x, ...) {
    show(x)
    invisible(x)
  }
)

