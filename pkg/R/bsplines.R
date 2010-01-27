bspline.basis <- function (x, nbasis, degree = 3) {
  .Call(bspline_basis,x,nbasis,degree)
}

periodic.bspline.basis <- function (x, nbasis, degree = 3, period = 1) {
  .Call(periodic_bspline_basis,x,nbasis,degree,period)
}
