##' accumulator variables
##'
##' Latent state variables that accumulate quantities through time.
##'
##' @name accumvars
##' @rdname accumvars
##' @family implementation information
##' @seealso \code{\link{sir}}
##' @details
##' In formulating models, one sometimes wishes to define a state variable that will accumulate some quantity over the interval between successive observations.
##' \pkg{pomp} provides a facility to make such features more convenient.
##' Specifically, if \eqn{a} is a state-variable named in the \code{pomp}'s \code{accumvars} argument, then for each interval \eqn{(t_k,t_{k+1})}{(t[k],t[k+1])}, \eqn{k=0,1,2,\dots}, \eqn{a} will be set to zero at prior to any \code{\link{rprocess}} computation over that interval.
##' Deterministic trajectory computation is handled slightly differently:
##' see \code{\link{flow}}.
##' For examples, see \code{\link{sir}} and the tutorials on the \href{https://kingaa.github.io/pomp/}{package website}.
##' @example examples/accumulators.R
##'
NULL
