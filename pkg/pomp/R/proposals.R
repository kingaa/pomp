mvn.rw.proposal.fn <- function (rw.sd) {
  parnm <- colnames(rw.sd)
  n <- ncol(rw.sd)
  if (is.null(parnm))
    stop(sQuote("rw.sd")," must have names")
  if (is.matrix(rw.sd)) {
    if (nrow(rw.sd)!=ncol(rw.sd))
      stop(sQuote("rw.sd")," must be a square matrix")
    ch <- try (chol(rw.sd,pivot=TRUE))
    if (inherits(ch,"try-error"))
      stop("error in Choleski factorization of ",sQuote("rw.sd"))
    oo <- order(attr(ch,"pivot"))
    Q <- ch[,oo]
    fn <- function (theta) {
      theta[parnm] <- theta[parnm]+Q.rnorm(n=n,mean=0,sd=1)
      theta
    }
  } else {
    fn <- function (theta) {
      theta[parnm] <- rnorm(n=n,mean=theta[parnm],sd=rw.sd)
      theta
    }
  }
  fn
}
