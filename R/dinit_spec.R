##' The initial-state distribution
##'
##' Specification of the initial-state distribution density evaluator, dinit.
##'
##' @name dinit specification
##' @rdname dinit_spec
##' @family implementation information
##' @seealso \code{\link{dinit}}
##' @details
##' To fully specify the unobserved Markov state process, one must give its distribution at the zero-time (\code{t0}).
##' One specifies how to evaluate the log probability density function for this distribution using the \code{dinit} argument.
##' As usual, this can be provided either as a C snippet or as an \R function.
##' In the former case, bear in mind that:
##' \enumerate{
##'   \item The goal of a this snippet is computation of a log likelihood, to be put into a variable named \code{loglik}.
##'   \item In addition to the state variables, parameters, and covariates (if any), the variable \code{t}, containing the zero-time, will be defined in the context in which the snippet is executed.
##' }
##' \link[=Csnippet]{General rules for writing C snippets can be found here}.
##'
##' If an \R function is to be used, pass
##' \preformatted{
##'    dinit = f
##' }
##' to \code{pomp}, where \code{f} is a function with arguments that can include the time \code{t}, any or all of the model state variables, parameters, and covariates.
##' As usual, \code{f} may take additional arguments, provided these are passed along with it in the call to \code{pomp}.
##' \code{f} must return a single numeric value, the log likelihood.
##'
##' @inheritSection pomp Note for Windows users
##'
NULL
