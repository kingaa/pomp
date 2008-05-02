reulermultinom <- function (n = 1, size, rate, dt) {
  ntrans <- length(rate)
  if (length(size)>1)
    warning("only the first element of 'size' is meansingful")
  size <- floor(size)
  x <- .C(
          "reulermultinom_multi",
          n = as.integer(n),
          ntrans=as.integer(ntrans),
          size=as.double(size),
          rate=as.double(rate),
          dt=as.double(dt),
          trans=double(n*ntrans),
          NAOK=FALSE,
          DUP=FALSE,
          PACKAGE="pomp"
          )$trans
  dim(x) <- c(ntrans,n)
  rownames(x) <- names(rate)
  x
}

deulermultinom <- function (x, size, rate, dt, log = FALSE) {
  ntrans <- length(rate)
  if (NROW(x)!=ntrans)
    stop("deulermultinom: length of 'x' should match length of 'rate'")
  n <- NCOL(x)
  if (length(size)>1)
    warning("only the first element of 'size' is meaningful")
  if (length(dt)>1)
    warning("only the first element of 'dt' is meaningful")
  size <- floor(size[1])
  f <- .C(
          "deulermultinom_multi",
          n=as.integer(n),
          ntrans=as.integer(ntrans),
          size=as.double(size),
          rate=as.double(rate),
          trans=as.double(x),
          dt=as.double(dt),
          give_log=as.integer(log),
          f=double(n),
          NAOK=FALSE,
          DUP=FALSE,
          PACKAGE="pomp"
          )$f
  f
}

