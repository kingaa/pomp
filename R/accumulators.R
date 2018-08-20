##' accumulators
##'
##' Accumulator variables
##'
##' @name accumulators
##' @rdname accumulators
##' @aliases zeronames
##' @family information on model implementation
##' @seealso \code{\link{sir}}
##'
##' @details
##' In formulating models, one sometimes wishes to define a state variable that will accumulate some quantity over the interval between successive observations.
##' \pkg{pomp} provides a facility to make such features more convenient.
##' Specifically, variables named in the \code{pomp}'s \code{zeronames} argument will be set to zero immediately following each observation.
##' See \code{\link{sir}} and the tutorials on the \href{https://kingaa.github.io/pomp/}{package website} for examples.
NULL
