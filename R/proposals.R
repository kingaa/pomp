##' MCMC proposal distributions
##'
##' Functions to construct proposal distributions for use with MCMC methods.
##'
##'
##' @name proposals
##' @rdname proposals
##'
##' @param rw.var square numeric matrix with row- and column-names.  Specifies
##' the variance-covariance matrix for a multivariate normal random-walk
##' proposal distribution.
##' @param rw.sd named numeric vector; random-walk SDs for a multivariate
##' normal random-walk proposal with diagonal variance-covariance matrix.
##' @param scale.start,scale.cooling,shape.start,target,max.scaling parameters
##' to control the proposal adaptation algorithm.  Beginning with MCMC
##' iteration \code{scale.start}, the scale of the proposal covariance matrix
##' will be adjusted in an effort to match the \code{target} acceptance ratio.
##' This initial scale adjustment is \dQuote{cooled}, i.e., the adjustment
##' diminishes as the chain moves along.  The parameter \code{scale.cooling}
##' specifies the cooling schedule: at n iterations after \code{scale.start},
##' the current scaling factor is multiplied with \code{scale.cooling^n}.  The
##' maximum scaling factor allowed at any one iteration is \code{max.scaling}.
##' After \code{shape.start} accepted proposals have accumulated, a scaled
##' empirical covariance matrix will be used for the proposals, following
##' Roberts and Rosenthal (2009).
##'
##' @return Each of these calls constructs a function suitable for use as the
##' \code{proposal} argument of \code{pmcmc} or \code{abc}.  Given a parameter
##' vector, each such function returns a single draw from the corresponding
##' proposal distribution.
##'
##' @author Aaron A. King, Sebastian Funk
##'
##' @seealso \code{\link{pmcmc}}, \code{\link{abc}}
##'
##' @references
##' Gareth O. Roberts and Jeffrey S. Rosenthal.
##' Examples of Adaptive MCMC.
##' J. Comput. Graph. Stat., 18:349--367, 2009.
NULL

##' @rdname proposals
mvn.diag.rw <- function (rw.sd) {
  if (missing(rw.sd) || !is.numeric(rw.sd)) {
    stop("in ",sQuote("mvn.diag.rw"),": ",
      sQuote("rw.sd")," must be a named numeric vector",call.=FALSE)
  }
  rw.sd <- rw.sd[rw.sd>0]
  parnm <- names(rw.sd)
  n <- length(rw.sd)
  if (is.null(parnm))
    stop("in ",sQuote("mvn.diag.rw"),": ",
      sQuote("rw.sd")," must have names",call.=FALSE)
  function (theta, ...) {
    theta[parnm] <- rnorm(n=n,mean=theta[parnm],sd=rw.sd)
    theta
  }
}

##' @rdname proposals
mvn.rw <- function (rw.var) {
  rw.var <- as.matrix(rw.var)
  parnm <- colnames(rw.var)
  ep <- paste0("in ",sQuote("mvn.rw"),": ")
  if (is.null(parnm))
    stop(ep,sQuote("rw.var")," must have row- and column-names",call.=FALSE)
  if (nrow(rw.var)!=ncol(rw.var))
    stop(ep,sQuote("rw.var")," must be a square matrix",call.=FALSE)
  ch <- chol(rw.var,pivot=TRUE)
  if (attr(ch,"rank") < ncol(ch))
    warning(ep,"rank-deficient covariance matrix",call.=FALSE)
  oo <- order(attr(ch,"pivot"))
  n <- Q <- NULL                 # to evade R CMD check false positive
  e <- new.env()
  e$Q <- ch[,oo]
  e$n <- ncol(rw.var)
  e$parnm <- parnm
  f <- function (theta, ...) {
    theta[parnm] <- theta[parnm]+rnorm(n=n,mean=0,sd=1)%*%Q
    theta
  }
  environment(f) <- e
  f
}

##' @rdname proposals
## a stateful function implementing an adaptive proposal
mvn.rw.adaptive <- function (rw.sd, rw.var, scale.start = NA,
  scale.cooling = 0.999,shape.start = NA, target = 0.234, max.scaling = 50) {

  ep <- paste0("in ",sQuote("mvn.rw.adaptive"),": ")

  if (!xor(missing(rw.sd),missing(rw.var))) {
    stop(ep,"exactly one of ",sQuote("rw.sd")," and ",sQuote("rw.var"),
      " must be supplied",call.=FALSE)
  }
  if (!missing(rw.var)) { ## variance supplied
    rw.var <- as.matrix(rw.var)
    parnm <- colnames(rw.var)
    if (is.null(parnm))
      stop(ep,sQuote("rw.var")," must have row- and column-names",call.=FALSE)
    if (nrow(rw.var)!=ncol(rw.var))
      stop(ep,sQuote("rw.var")," must be a square matrix",call.=FALSE)
    if (any(parnm!=rownames(rw.var)))
      stop(ep,"row- and column-names of ",sQuote("rw.var"),
        " must agree",call.=FALSE)
  } else if (!missing(rw.sd)) { ## sd supplied (diagonal)
    if (!is.numeric(rw.sd) || is.null(names(rw.sd))) {
      stop(ep,sQuote("rw.sd")," must be a named numeric vector",call.=FALSE)
    }
    rw.sd <- rw.sd[rw.sd>0]
    parnm <- names(rw.sd)
    rw.var <- diag(rw.sd^2,nrow=length(rw.sd),ncol=length(rw.sd))
    dimnames(rw.var) <- list(parnm,parnm)
  }

  scale.start <- as.integer(scale.start)
  scale.cooling <- as.numeric(scale.cooling)
  shape.start <- as.integer(shape.start)
  target <- as.numeric(target) ## target acceptance ratio
  if (!isTRUE(scale.start > 0))
    stop(ep,sQuote("scale.start")," must be a positive integer",call.=FALSE)
  if (!isTRUE((scale.cooling > 0) && (scale.cooling <= 1)))
    stop(ep,sQuote("scale.cooling")," must be in (0,1]",call.=FALSE)
  if (!isTRUE(shape.start > 0))
    stop(ep,sQuote("shape.start")," must be a positive integer",call.=FALSE)
  if (!isTRUE(target > 0 && target < 1)) {
    stop(ep,sQuote("target")," must be a number in (0,1)",call.=FALSE)
  }

  ## variables that will follow 'f'
  scaling <- 1
  theta.mean <- NULL
  covmat.emp <- array(data=0,dim=dim(rw.var),dimnames=dimnames(rw.var))

  function (theta, .n, .accepts, verbose, ...) {
    if (.n == 0) return(theta) ## handle initial test run by pmcmc
    if (is.null(theta.mean)) theta.mean <<- theta[parnm]
    if (!is.na(scale.start) && .n >= scale.start &&
        (is.na(shape.start) || .accepts < shape.start)) {
      ## adapt size of covmat until we get enough accepted jumps
      scaling <<- min(scaling*exp(scale.cooling^(.n-scale.start)*
          (.accepts/.n-target)),
        max.scaling)
      covmat <- scaling^2*rw.var
    } else if (!is.na(shape.start) && .accepts >= shape.start) {
      scaling <<- 2.38^2/length(parnm)
      covmat <- scaling*covmat.emp
    } else {
      covmat <- rw.var
    }
    if (verbose) {
      cat("proposal covariance matrix:\n")
      print(covmat)
    }
    if (!is.na(shape.start)) {
      theta.mean <<- ((.n-1)*theta.mean+theta[parnm])/.n
      covmat.emp <<- ((.n-1)*covmat.emp+tcrossprod(theta[parnm]-theta.mean))/.n
    }
    ch <- chol(covmat,pivot=TRUE)
    if (attr(ch,"rank")<length(parnm))
      warning(ep,"degenerate proposal",call.=FALSE)
    oo <- order(attr(ch,"pivot"))
    Q <- ch[,oo]
    theta[parnm] <- theta[parnm]+rnorm(n=length(parnm),mean=0,sd=1)%*%Q
    theta
  }
}
