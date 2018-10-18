##' @include listie.R
##' @keywords internal
##' @rdname listie
#' @name [-listie
#' @aliases [,listie-method
NULL

##' @export
##' @rdname listie
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
      new(cl,y)
    }
  }
)
