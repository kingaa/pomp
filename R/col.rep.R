col.rep <- function (x, n, named = TRUE) {
  if(is.vector(x)) {
    y <- matrix(x,nrow=length(x),ncol=n)
    if (named && is.null(names(x))) stop("x should be a named vector")
    if (named) rownames(y) <- names(x)
  }
  if (is.matrix(x)) {
    if(dim(x)[2]>1) stop("x should be an m x 1 matrix")
    y <- matrix(x,nrow=dim(x)[1],ncol=n)
    if (named && is.null(rownames(x))) stop("x should have rownames")
    if (named) rownames(y) <- rownames(x)
  }
  y
}


