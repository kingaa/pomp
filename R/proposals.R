mvn.diag.rw <- function (rw.sd) {
  if (!is.numeric(rw.sd)) {
    stop(sQuote("rw.sd")," must be a named numeric vector")
  }
  rw.sd <- rw.sd[rw.sd>0]
  parnm <- names(rw.sd)
  n <- length(rw.sd)
  if (is.null(parnm))
    stop(sQuote("rw.sd")," must have names")
  function (theta) {
    theta[parnm] <- rnorm(n=n,mean=theta[parnm],sd=rw.sd)
    theta
  }
}

mvn.rw <- function (rw.var) {
  rw.var <- as.matrix(rw.var)
  parnm <- colnames(rw.var)
  n <- ncol(rw.var)
  if (is.null(parnm))
    stop(sQuote("rw.var")," must have row- and column-names")
  if (nrow(rw.var)!=ncol(rw.var))
    stop(sQuote("rw.var")," must be a square matrix")
  ch <- try (chol(rw.var,pivot=TRUE))
  if (inherits(ch,"try-error"))
    stop("error in Choleski factorization of ",sQuote("rw.var"))
  oo <- order(attr(ch,"pivot"))
  Q <- ch[,oo]
  function (theta) {
    theta[parnm] <- theta[parnm]+rnorm(n=n,mean=0,sd=1)%*%Q
    theta
  }
}
