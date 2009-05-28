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
  if (nbasis<degree)
    stop("periodic.bspline.basis error: must have ",sQuote("nbasis")," >= ",sQuote("degree"),call.=FALSE)
  dx <- period/nbasis
  knots <- seq(-degree*dx,period+degree*dx,by=dx)
  y <- .Call(bspline_basis,x%%period,degree,knots)
  if (degree>0)
    y[,1:degree] <- y[,1:degree]+y[,-(1:nbasis)]
  shift <- floor((degree-1)/2)
  y[,c(seq(from=shift+1,to=nbasis,by=1),seq_len(shift))]
}
