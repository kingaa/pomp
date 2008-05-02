bspline <- function (x, i, degree, knots)
  .Call(bspline_basis_function,x,as.integer(i),as.integer(degree),knots)

bspline.basis <- function (x, degree = 3, knots)
  .Call(bspline_basis,x,as.integer(degree),knots)

periodic.bspline.basis <- function (x, nbasis, degree = 3, period = 1) {
  if (nbasis<degree)
    stop("must have nbasis >= degree")
  dx <- period/nbasis
  knots <- seq(-degree*dx,period+degree*dx,by=dx)
  y <- bspline.basis(x%%period,degree,knots)
  if (degree>0)
    y[,1:degree] <- y[,1:degree]+y[,-(1:nbasis)]
  shift <- floor((degree-1)/2)
  y[,c(seq(from=shift+1,to=nbasis,by=1),seq_len(shift))]
}
