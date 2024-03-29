##' dmeasure specification
##'
##' Specification of the measurement model density function, dmeasure.
##'
##' @name dmeasure_spec
##' @rdname dmeasure_spec
##' @family implementation information
##' @seealso \code{\link{dmeasure}}
##' @details
##' The measurement model is the link between the data and the unobserved state process.
##' It can be specified either by using one or both of the \code{rmeasure} and \code{dmeasure} arguments.
##'
##' Suppose you have a procedure to compute the probability density of an observation given the value of the latent state variables.
##' Then you can furnish \preformatted{
##'    dmeasure = f}
##' to \pkg{pomp} algorithms,
##' where \code{f} is a C snippet or \R function that implements your procedure.
##'
##' Using a C snippet is much preferred, due to its much greater computational efficiency.
##' See \code{\link{Csnippet}} for general rules on writing C snippets.
##' The goal of a \dfn{dmeasure} C snippet is to fill the variable \code{lik} with the either the probability density or the log probability density, depending on the value of the variable \code{give_log}.
##'
##' In writing a \code{dmeasure} C snippet, observe that:
##'   \enumerate{
##'     \item In addition to the states, parameters, covariates (if any), and observables, the variable \code{t}, containing the time of the observation will be defined in the context in which the snippet is executed.
##'     \item Moreover, the Boolean variable \code{give_log} will be defined.
##'     \item The goal of a dmeasure C snippet is to set the value of the \code{lik} variable to the likelihood of the data given the state, if \code{give_log == 0}.
##'     If \code{give_log == 1}, \code{lik} should be set to the log likelihood.
##'   }
##'
##' If \code{dmeasure} is to be provided instead as an \R function, this is accomplished by supplying \preformatted{
##'   dmeasure = f}
##' to \code{pomp}, where \code{f} is a function.
##' The arguments of \code{f} should be chosen from among the observables, state variables, parameters, covariates, and time.
##' It must also have the arguments \code{\dots}, and \code{log}.
##' It can take additional arguments via the \link[=userdata]{userdata facility}.
##' \code{f} must return a single numeric value, the probability density (or log probability density if \code{log = TRUE}) of \code{y} given \code{x} at time \code{t}.
##'
##' @example examples/dmeasure_spec.R
##'
##' @section Important note:
##' \strong{It is a common error to fail to account for both \code{log = TRUE} and \code{log = FALSE} when writing the \code{dmeasure} C snippet or function.}
##'
##' @section Default behavior:
##' If \code{dmeasure} is left unspecified, calls to \code{\link{dmeasure}} will return missing values (\code{NA}).
##'
##' @inheritSection pomp Note for Windows users
##'
NULL
