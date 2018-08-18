##' @include listies.R
##' @keywords internal
##' @rdname listies
#' @name [-listies
#' @aliases [,listies-method
NULL

##' @export
##' @rdname listies
setMethod(
  "[",
  signature=signature(x="listies"),
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
