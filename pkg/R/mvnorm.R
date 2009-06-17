## multivariate normal density
## stolen from package 'mvtnorm' version 0.9-7


rmvnorm <- function (n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)),
                     method=c("eigen", "svd", "chol"))
{    
    if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps))) {
        stop("sigma must be a symmetric matrix")
    }
    if (length(mean) != nrow(sigma)) {
        stop("mean and sigma have non-conforming size")
    }
    sigma1 <- sigma
    dimnames(sigma1) <- NULL
    if(!isTRUE(all.equal(sigma1, t(sigma1)))){
        warning("sigma is numerically not symmetric")
    }

    method <- match.arg(method)
    
    if(method == "eigen"){
        ev <- eigen(sigma, symmetric = TRUE)
        if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))){
            warning("sigma is numerically not positive definite")
        }    
        retval <- ev$vectors %*%  diag(sqrt(ev$values),length(ev$values)) %*% t(ev$vectors)
    }
    else if(method == "svd"){
        sigsvd <- svd(sigma)
        if (!all(sigsvd$d >= -sqrt(.Machine$double.eps) * abs(sigsvd$d[1]))){
            warning("sigma is numerically not positive definite")
        }    
        retval <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
    }    
    else if(method == "chol"){
        retval <- chol(sigma, pivot = TRUE)
        o <- order(attr(retval, "pivot"))
        retval <- retval[,o]
    }
    
    retval <- matrix(rnorm(n * ncol(sigma)), nrow = n) %*%  retval
    retval <- sweep(retval, 2, mean, "+")
    retval
}

dmvnorm <- function (x, mean, sigma, log=FALSE)
{
    if (is.vector(x)) {
        x <- matrix(x, ncol = length(x))
    }
    if (missing(mean)) {
        mean <- rep(0, length = ncol(x))
    }
    if (missing(sigma)) {
        sigma <- diag(ncol(x))
    }
    if (NCOL(x) != NCOL(sigma)) {
        stop("x and sigma have non-conforming size")
    }
    if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps))) {
        stop("sigma must be a symmetric matrix")
    }
    if (length(mean) != NROW(sigma)) {
        stop("mean and sigma have non-conforming size")
    }
    distval <- mahalanobis(x, center = mean, cov = sigma)
    logdet <- sum(log(eigen(sigma, symmetric=TRUE, only.values=TRUE)$values))
    logretval <- -(ncol(x)*log(2*pi) + logdet + distval)/2
    if(log) return(logretval)
    exp(logretval)
}
