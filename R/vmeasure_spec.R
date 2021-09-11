##' The variance of the measurement model
##'
##' Specification of the measurement-model covariance  matrix, vmeasure.
##'
##' @name vmeasure specification
##' @rdname vmeasure_spec
##' @family implementation information
##' @seealso \code{\link{vmeasure}}
##' 
##' @details
##' The measurement model is the link between the data and the unobserved state process.
##' Some algorithms require the conditional covariance of the measurement model, given the latent state and parameters.
##' This is supplied using the \code{vmeasure} argument.
##'
##' Suppose you have a procedure to compute this conditional covariance matrix, given the value of the latent state variables.
##' Then you can furnish \preformatted{
##'   vmeasure = f}
##' to \pkg{pomp} algorithms,
##' where \code{f} is a C snippet or \R function that implements your procedure.
##'
##' Using a C snippet is much preferred, due to its much greater computational efficiency.
##' See \code{\link{Csnippet}} for general rules on writing C snippets.
##'
##' In writing a \code{vmeasure} C snippet, bear in mind that:
##'   \enumerate{
##'     \item The goal of such a snippet is to fill variables named \code{V_y_z} with the conditional covariances of observables \code{y}, \code{z}.
##'     Accordingly, there should be one assignment of \code{V_y_z} and one assignment of \code{V_z_y} for each pair of observables \code{y} and \code{z}.
##'     \item In addition to the states, parameters, and covariates (if any), the variable \code{t}, containing the time of the observation, will be defined in the context in which the snippet is executed.
##'   }
##'
##' The demos and the tutorials on the \href{https://kingaa.github.io/pomp/}{package website} give examples.
##'
##' It is also possible, though less efficient, to specify \code{vmeasure} using an \R function.
##' In this case, specify it by furnishing \preformatted{
##'   vmeasure = f}
##' to \code{pomp}, where \code{f} is an \R function.
##' The arguments of \code{f} should be chosen from among the state variables, parameters, covariates, and time.
##' It must also have the argument \code{...}.
##' \code{f} must return a square matrix of dimension equal to the number of observable variables.
##' The row- and column-names of this matrix should match the names of the observable variables.
##' The matrix should of course be symmetric.
##'
##' @section Default behavior:
##' The default \code{vmeasure} is undefined.
##' It will yield missing values (\code{NA}).
##' 
##' @inheritSection pomp Note for Windows users
##'
NULL
