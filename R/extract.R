##' @include listie.R
##' @keywords internal
##' @rdname listie
##' @name [-listie
##' @aliases [,listie-method
NULL

##' @rdname listie
##' @importFrom stats setNames
##' @export
setMethod(
  "[",
  signature=signature(x="listie"),
  definition=function(x, i, ...) {
    y <- as(x,"list")
    names(y) <- names(x)
    cl <- class(x)
    y <- unlist(y[i])
    if (is.null(y)) {
      list(NULL)
    } else {
      setNames(new(cl,y),names(y))
    }
  }
)
