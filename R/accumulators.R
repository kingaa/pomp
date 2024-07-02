##' accumulator variables
##'
##' Latent state variables that accumulate quantities through time.
##'
##' @name accumvars
##' @rdname accumvars
##' @family implementation information
##' @seealso \code{\link{sir}}
##'
##' @details
##' In formulating models, one sometimes wishes to define a state variable that will accumulate some quantity over the interval between successive observations.
##' \pkg{pomp} provides a facility to make such features more convenient.
##' Specifically, variables named in the \code{pomp}'s \code{accumvars} argument will be set to zero immediately following each observation.
##' See \code{\link{sir}} and the tutorials on the \href{https://kingaa.github.io/pomp/}{package website} for examples.
##'
##' @example examples/accumulators.R
##'
NULL
