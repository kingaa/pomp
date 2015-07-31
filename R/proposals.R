mvn.diag.rw <- function (rw.sd, log = FALSE) {
  if (!is.numeric(rw.sd)) {
    stop(sQuote("rw.sd")," must be a named numeric vector")
  }
  rw.sd <- rw.sd[rw.sd>0]
  parnm <- names(rw.sd)
  n <- length(rw.sd)
  if (is.null(parnm))
    stop(sQuote("rw.sd")," must have names")
  if (log) {
    function (theta) {
      theta[parnm] <- rlnorm(n=n,meanlog=log(theta[parnm]),
                             sdlog=rw.sd)
      theta
    }
  } else {
    function (theta) {
      theta[parnm] <- rnorm(n=n,mean=theta[parnm],sd=rw.sd)
      theta
    }
  }
}

mvn.rw <- function (rw.var, log = FALSE) {
  log <- as.logical(log)
  rw.var <- as.matrix(rw.var)
  parnm <- colnames(rw.var)
  if (is.null(parnm))
    stop(sQuote("rw.var")," must have row- and column-names")
  if (nrow(rw.var)!=ncol(rw.var))
    stop(sQuote("rw.var")," must be a square matrix")
  ch <- try (chol(rw.var,pivot=TRUE))
  if (inherits(ch,"try-error"))
    stop("in ",sQuote("mvn.rw"),": error in Choleski factorization",
         call.=FALSE)
  if (attr(ch,"rank") < ncol(ch))
    warning("in ",sQuote("mvn.rw"),": rank-deficient covariance matrix",
            call.=FALSE)
  oo <- order(attr(ch,"pivot"))
  n <- Q <- NULL                 # to evade R CMD check false positive
  e <- new.env()
  e$Q <- ch[,oo]
  e$n <- ncol(rw.var)
  e$parnm <- parnm
  if (log) {
    f <- function (theta) {
      theta[parnm] <- theta[parnm]*exp(rnorm(n=n,mean=0,sd=1)%*%Q)
      theta
    }
  } else {
    f <- function (theta) {
      theta[parnm] <- theta[parnm]+rnorm(n=n,mean=0,sd=1)%*%Q
      theta
    }
  }
  environment(f) <- e
  f
}
