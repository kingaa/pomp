bspline <- function (x, i, degree, knots)
  .Call(bspline_basis_function,x,as.integer(i),as.integer(degree),knots)

bspline.basis <- function (x, degree = 3, knots)
  .Call(bspline_basis,x,as.integer(degree),knots)

periodic.bspline.basis <- function (x, nbasis, degree = 3, period = 1) {
  if (any(x < 0) || any(x > period))
    stop("cannot evaluate the basis outside the fundamental domain")
  if (nbasis < degree)
    stop("must have nbasis >= degree")
  dx <- period/nbasis
  knots <- seq(-degree*dx,period+degree*dx,by=dx)
  y <- bspline.basis(x,degree,knots)
  if (degree > 0)
    y[,1:degree] <- y[,1:degree]+y[,-(1:nbasis)]
  shift <- floor((degree-1)/2)
  y[,c((shift+1):nbasis,1:shift)]
}
