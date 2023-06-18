##' Beta-binomial distribution
##'
##' Density and random generation for the Beta-binomial distribution with parameters \code{size}, \code{mu}, and \code{theta}.
##' 
##' A variable \eqn{X} is Beta-binomially distributed if
##' \eqn{X\sim{\mathrm{Binomial}(n,P)}}{X~{Binomial(n,P)}} where \eqn{P\sim{\mathrm{Beta}(\mu,\theta)}}{P~Beta(mu,theta)}.
##' Using the standard \eqn{(a,b)} parameterization, \eqn{a=\mu\,\theta}{a=mu*theta} and \eqn{b=(1-\mu)\,\theta}{b=(1-mu)*theta}.
##'
##' @name betabinomial
##' @rdname betabinom
##' @family implementation information
##' @concept probability distributions
##' @inheritSection eulermultinom C API
##' @inheritParams eulermultinom
##' @param size \code{size} parameter of the binomial distribution
##' @param prob mean of the Beta distribution
##' @param theta Beta distribution dispersion parameter
##' @param x vector of non-negative integer quantiles
##' @return
##' \item{rbetabinom}{
##'    Returns a vector of length \code{n} containing random variates drawn from the Beta-binomial distribution.
##' }
##' \item{dbetabinom}{
##'    Returns a vector (of length equal to the number of columns of \code{x}) containing the probabilities of observing each column of \code{x} given the specified parameters (\code{size}, \code{prob}, \code{theta}).
##' }
NULL

##' @rdname betabinom
##' @export
rbetabinom <- function (n = 1, size, prob, theta) {
  tryCatch(
    .Call(P_R_BetaBinom,n,size,prob,theta),
    error = function (e) pStop(who="rbetabinom",conditionMessage(e))
  )
}

##' @rdname betabinom
##' @export
dbetabinom <- function (x, size, prob, theta, log = FALSE) {
  tryCatch(
    .Call(P_D_BetaBinom,x,size,prob,theta,log),
    error = function (e) pStop(who="dbetabinom",conditionMessage(e))
  )
}
