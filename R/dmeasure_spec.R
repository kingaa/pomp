##' The measurement model density
##'
##' Specification of dmeasure.
##'
##' @name dmeasure_spec
##' @rdname dmeasure_spec
##' @family information on model implementation
##'
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
##' to \code{pomp}, where \code{f} has prototype \preformatted{
##'   f(y, x, t, params, log, \dots)}
##' Again, it can take additional arguments that are passed with it in the call to \code{pomp}.
##' When \code{f} is called,
##' \itemize{
##'   \item \code{y} will be a named numeric vector of length \code{nobs} containing values of the observed variables;
##'   \item \code{x} will be a named numeric vector of length \code{nvar} containing state variables;
##'   \item \code{params} will be a named numeric vector of length \code{npar} containing parameters;
##'   \item \code{t} will be a scalar, the corresponding observation time.
##' }
##' \code{f} must return a single numeric value, the probability density (or log probability density) of \code{y} given \code{x} at time \code{t}.
##'
##' \strong{NB: it is a common error to fail to account for both \code{log = TRUE} and \code{log = FALSE} when writing the \code{dmeasure} C snippet or function.}
##'
NULL