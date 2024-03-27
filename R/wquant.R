##' Weighted quantile function
##'
##' Estimate weighted quantiles.
##'
##' \code{wquant} estimates quantiles of weighted data using the estimator of Harrell & Davis (1982), with improvements recommended by Andrey Akinshin.
##'
##' @param x numeric; a vector of data.
##' @param weights numeric; vector of weights.
##' @param probs numeric; desired quantiles.
##'
##' @return
##' \code{wquant} returns a vector containing the estimated quantiles.
##' If \code{probs} has names, these are inherited.
##'
##' @author Aaron A. King
##' @include package.R
##' @examples
##' x <- c(1,1,1,2,2,3,3,3,3,4,5,5,6,6,6)
##' quantile(x)
##' wquant(x)
##' wquant(c(1,2,3,4,5,6),weights=c(3,2,4,1,2,3))
##' wquant(c(1,2,3,4,5),c(1,0,0,1,1))
##' @references
##' \Harrell1982
##'
## The discussion in Andrey Akinshin's blog post
## https://aakinshin.net/posts/weighted-quantiles/
## was instrumental in developing these codes.
##'
##' @rdname wquant
##' @importFrom stats pbeta
##' @export
wquant <- function (
  x, weights = rep(1, length(x)),
  probs = c(`0%`=0, `25%`=0.25, `50%`=0.5, `75%`=0.75, `100%`=1)
) {
  x <- as.numeric(x)
  weights <- as.numeric(weights)
  if (length(weights) != length(x))
    pStop(sQuote("x")," and ",sQuote("weights")," must be of equal length.")
  if (any(is.na(x)) || any(!is.finite(weights)))
    pStop("NA and non-finite values are disallowed.")
  if (any(weights < 0))
    pStop("weights must be non-negative.")
  if (!is.numeric(probs) || any(is.na(probs)) ||
        isTRUE(any(probs < 0 | probs > 1))) {
    pStop(sQuote("probs")," must be a numeric vector with values in [0,1].")
  }
  ## order the data and the weights
  idx <- order(x)
  x <- x[idx]
  weights <- weights[idx]
  ## accumulate and normalize the weights
  w <- cumsum(c(0,weights))
  w <- w/w[length(w)]
  q <- probs       # inherits the names of `probs`
  ess <- sum(weights)^2/sum(weights^2) # Kish effective sample size
  a <- probs*(ess+1) # a,b are shape parameters for Beta distribution
  b <- (1-probs)*(ess+1)
  for (j in seq_along(probs)) {
    W <- pbeta(q=w,shape1=a[j],shape2=b[j])
    W[w==0] <- 0
    W[w==1] <- 1
    q[j] <- sum(diff(W)*x)
  }
  q
}
