##' The expectation of the measurement model
##'
##' Specification of the measurement-model conditional expectation, emeasure.
##'
##' @name emeasure specification
##' @rdname emeasure_spec
##' @family implementation information
##' @seealso \code{\link{emeasure}}
##' 
##' @details
##' The measurement model is the link between the data and the unobserved state process.
##' Some algorithms require the conditional expectation of the measurement model, given the latent state and parameters.
##' This is supplied using the \code{emeasure} argument.
##'
##' Suppose you have a procedure to compute this conditional expectation, given the value of the latent state variables.
##' Then you can furnish \preformatted{
##'   emeasure = f}
##' to \pkg{pomp} algorithms,
##' where \code{f} is a C snippet or \R function that implements your procedure.
##'
##' Using a C snippet is much preferred, due to its much greater computational efficiency.
##' See \code{\link{Csnippet}} for general rules on writing C snippets.
##'
##' In writing an \code{emeasure} C snippet, bear in mind that:
##'   \enumerate{
##'     \item The goal of such a snippet is to fill variables named \code{E_y} with the conditional expectations of observables \code{y}.
##'     Accordingly, there should be one assignment of \code{E_y} for each observable \code{y}.
##'     \item In addition to the states, parameters, and covariates (if any), the variable \code{t}, containing the time of the observation, will be defined in the context in which the snippet is executed.
##'   }
##'
##' The demos and the tutorials on the \href{https://kingaa.github.io/pomp/}{package website} give examples.
##'
##' It is also possible, though less efficient, to specify \code{emeasure} using an \R function.
##' In this case, specify the measurement model expectation by furnishing \preformatted{
##'   emeasure = f}
##' to \code{pomp}, where \code{f} is an \R function.
##' The arguments of \code{f} should be chosen from among the state variables, parameters, covariates, and time.
##' It must also have the argument \code{...}.
##' \code{f} must return a named numeric vector of length equal to the number of observable variables.
##' The names should match those of the observable variables.
##'
##' @section Default behavior:
##' The default \code{emeasure} is undefined.
##' It will yield missing values (\code{NA}).
##' 
##' @inheritSection pomp Note for Windows users
##'
NULL
