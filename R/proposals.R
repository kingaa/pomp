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

## a stateful function implementing an adaptive proposal
mvn.rw.adaptive <- function (rw.sd, rw.var,
                             memory = 100, target = 0.234,
                             ...) {
  if (!xor(missing(rw.sd),missing(rw.var))) {
    stop(sQuote("mvn.rw.adaptive")," error: exactly one of ",
         sQuote("rw.sd")," and ",sQuote("rw.var"),
         " must be supplied",call.=FALSE)
  }
  if (!missing(rw.var)) { ## variance supplied
    rw.var <- as.matrix(rw.var)
    parnm <- colnames(rw.var)
    if (is.null(parnm))
      stop(sQuote("rw.var")," must have row- and column-names")
    if (any(parnm!=rownames(rw.var)))
      stop("row- and column-names of ",sQuote("rw.var"),
           " must agree")
    if (nrow(rw.var)!=ncol(rw.var))
      stop(sQuote("rw.var")," must be a square matrix")
  } else if (!missing(rw.sd)) { ## sd supplied (diagonal)
    if (!is.numeric(rw.sd) || is.null(names(rw.sd))) {
      stop(sQuote("rw.sd")," must be a named numeric vector")
    }
    rw.sd <- rw.sd[rw.sd>0]
    parnm <- names(rw.sd)
    rw.var <- diag(rw.sd^2)
    dimnames(rw.var) <- list(parnm,parnm)
  }

  memory <- as.integer(memory) ## duration of memory
  if (memory < 0) {
    stop(sQuote("mvn.rw.adaptive")," error: ",
         sQuote("memory")," must be a positive integer",call.=FALSE)
  }
  target <- as.numeric(target) ## target acceptance ratio
  if (target <= 0 || target >= 1) {
    stop(sQuote("mvn.rw.adaptive")," error: ",
         sQuote("target")," must be a number in (0,1)",call.=FALSE)
  }
  
  ## create environment to hold state of proposal distribution
  e <- new.env()
  e$covmat <- rw.var
  e$memory <- memory
  e$parnm <- parnm
  e$target <- target
  e$scaling <- 2
  e$theta.mean <- NULL
  e$acc.ratio <- NA
  e$naccept <- NA
  
  f <- function (theta, .n, .traces, .accepts, ...) {
    if (.n == 0) return(theta) ## to handle initial test that pmcmc runs
    if (.n >= memory) {    ## do adaptation
      if (.n == memory) {  ## look at empirical mean and variance and acceptance ratio
        assign(theta.mean <<- colMeans(.traces[1:memory,parnm])
        covmat <<- var(.traces[1:memory,parnm])
        acc.ratio <<- .accepts/.n
      } else { ## update mean, variance, and acceptance ratio
        theta.mean <<- ((memory-1)*theta.mean+theta[parnm])/memory
        covmat <<- ((memory-1)*covmat+tcrossprod(theta[parnm]-theta.mean))/memory
        acc.ratio <<- ((memory-1)*acc.ratio+(.accepts>naccept))/memory
      }
      naccept <<- .accepts # store this so we can see when we have an acceptance
      scaling <<- min(exp(acc.ratio-target),50) # adjustment to size of proposals
      covmat <<- scaling^2*covmat
    }
    if (verbose) {
      cat("proposal covariance matrix:\n")
      print(covmat)
    }
    ch <- try (chol(covmat,pivot=TRUE))
    ## we should worry more about degeneracy than this:
    if (inherits(ch,"try-error")) stop("Choleski factorization problem")
    if (attr(ch,"rank")<length(parnm)) stop("degenerate proposal")
    oo <- order(attr(ch,"pivot"))
    Q <- ch[,oo]
    theta[parnm] <- theta[parnm]+rnorm(n=length(parnm),mean=0,sd=1)%*%Q
    theta
  }

  environment(f) <- e
  f
}
