bspline <- function (x, i, degree, knots)
  .Call(bspline_basis_function,x,as.integer(i),as.integer(degree),knots)

bspline.basis <- function (x, nbasis, degree = 3) {
  if (nbasis<=degree)
    stop("bspline.basis error: must have ",sQuote("nbasis")," > ",sQuote("degree"),call.=FALSE)
  min.x <- min(x,na.rm=TRUE)
  max.x <- max(x,na.rm=TRUE)
  dx <- (max.x-min.x)/(nbasis-degree)
  tails <- degree*dx
  knots <- seq(from=min.x-tails,to=max.x+tails,by=dx)
  .Call(bspline_basis,x,degree,knots)
}

periodic.bspline.basis <- function (x, nbasis, degree = 3, period = 1) {
  .Call(periodic_bspline_basis,x,nbasis,degree,period)
}
