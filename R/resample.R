##' Resample
##'
##' Systematic resampling.
##'
##' @return
##' A vector of integers containing the indices of the resample.
##'
##' @rdname resample
##' @name resample
##' @keywords internal
##' @concept sampling
##'
##' @param weights numeric; vector of weights.
##' @param Np integer scalar; number of samples to draw.
##'
NULL

##' @export
##' @name systematic_resample
##' @rdname resample
##'
systematic_resample <- function (weights, Np = length(weights))
  .Call(P_systematic_resampling,weights,Np)
