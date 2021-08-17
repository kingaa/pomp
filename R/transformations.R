##' Transformations
##'
##' Some useful parameter transformations.
##'
##' Parameter transformations can be used in many cases to recast constrained optimization problems as unconstrained problems.
##' Although there are no limits to the transformations one can implement using the \code{\link{parameter_trans}} facilty, \pkg{pomp} provides a few ready-built functions to implement some very commonly useful ones.
##'
##' @name transformations
##' @rdname transformations
##' @family implementation information
##' @concept parameter transformations
##'
NULL

##' @export
##' @rdname transformations
##'
##' @param p numeric; a quantity in [0,1].
##'
##' @details
##' The logit transformation takes a probability \eqn{p} to its log odds, \eqn{\log\frac{p}{1-p}}{log(p/(1-p))}.
##' It maps the unit interval \eqn{[0,1]} into the extended real line \eqn{[-\infty,\infty]}.
##'
logit <- function (p) {
  .Call(P_LogitTransform,p)
}

##' @export
##' @rdname transformations
##'
##' @param x numeric; the log odds ratio.
##'
##' @details The inverse of the logit transformation is the expit transformation.
##'
expit <- function (x) {
  .Call(P_ExpitTransform,x)
}

##' @export
##' @rdname transformations
##'
##' @param X numeric; a vector containing the quantities to be transformed according to the log-barycentric transformation.
##'
##' @details
##' The log-barycentric transformation takes a vector \eqn{X_i}{Xi}, \eqn{i=1,\dots,n}, to a vector \eqn{Y_i}{Yi}, where \deqn{Y_i = \log\frac{X_i}{\sum_j X_j}.}{Yi = log(Xi/sum(Xj)).}
##' If \eqn{X} is an \eqn{n}-vector, it takes every simplex defined by \eqn{\sum_i X_i = c}{sum(Xi)=c}, \eqn{c} constant, to n-dimensional Euclidean space \eqn{R^n}{R^n}.
##'
log_barycentric <- function (X) {
  .Call(P_LogBarycentricTransform,X)
}

##' @export
##' @rdname transformations
##'
##' @param Y numeric; a vector containing the log fractions.
##'
##' @details
##' The inverse of the log-barycentric transformation is implemented as \code{inv_log_barycentric}.
##' Note that it is not a true inverse, in the sense that it takes \eqn{R^n} to the \emph{unit} simplex, \eqn{\sum_i X_i = 1}{sum(Xi)=1}.
##' Thus, \preformatted{
##'     log_barycentric(inv_log_barycentric(Y)) == Y,
##' } but \preformatted{
##'     inv_log_barycentric(log_barycentric(X)) == X
##' } only if \code{sum(X) == 1}.
##'
inv_log_barycentric <- function (Y) {
  .Call(P_InverseLogBarycentricTransform,Y)
}
