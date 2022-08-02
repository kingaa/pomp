##' Weighted quantile function
##'
##' Computes weighted quantiles.
##'
##' \code{wquant} computes a quantile of type 7 according to the typology of \code{\link[stats]{quantile}}.
##' 
##' @param x numeric; a vector of data.
##' @param weights numeric; vector of weights.
##' @param probs numeric; desired quantiles.
##'
##' @return
##' \code{wquant} returns a vector containing the quantiles.
##' If \code{probs} has names, these are inherited.
##' 
##' @author Aaron A. King
##' @include package.R
##' @examples
##' x <- c(1,1,1,2,2,3,3,3,3,4,5,5,6,6,6)
##' quantile(x)
##' wquant(c(1,2,3,4,5,6),weights=c(3,2,4,1,2,3))
##'
##' \dontshow{
##'   stopifnot(quantile(x)==wquant(c(1,2,3,4,5,6),weights=c(3,2,4,1,2,3)))
##'   try(wquant(c(1,NA),c(1,2)))
##'   try(wquant(c(1,2),c(NA,1)))
##'   try(wquant(c(1,2,3),c(1,2)))
##'   try(wquant(c(1,2,3),c(1,1,1),probs=c(0.1,NA)))
##'   try(wquant(c(1,2,3),c(1,2,3),probs=c(0.1,2)))
##' }
##'
##' @rdname wquant
##' @importFrom stats approx setNames
##' @export
wquant <- function (
  x, weights,
  probs = c(`0%`=0, `25%`=0.25, `50%`=0.5, `75%`=0.75, `100%`=1)
) {
  x <- as.numeric(x)
  weights <- as.numeric(weights)
  if (length(x)!=length(weights))
    pStop("wquant",sQuote("x")," and ",sQuote("weights"),
      " must be of equal length.")
  if (any(is.na(x)) || any(is.na(weights)))
    pStop("wquant","NA values are disallowed.")
  if (!is.numeric(probs) || any(is.na(probs)) ||
        isTRUE(any(probs < 0 | probs > 1))) {
    pStop("wquant",sQuote("probs"),
      " must be a numeric vector with values in [0,1].")
  }
  idx <- order(x)
  x <- x[idx]
  weights <- weights[idx]
  n <- sum(weights)
  k <- length(probs)
  ord <- 1+(n-1)*probs
  low <- pmax(floor(ord),1)
  high <- pmin(low+1,n)
  ord <- ord%%1
  allq <- approx(
    x=cumsum(weights),y=x,xout=c(low,high),
    method="constant",f=1,rule=2
  )$y
  dim(allq) <- c(k,2)
  qs <- (1-ord)*allq[,1]+ord*allq[,2]
  setNames(qs,names(probs))
}
