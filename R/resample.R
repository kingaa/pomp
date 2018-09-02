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
##'
NULL

##' @export
##' @name systematic_resample
##' @rdname resample
##'
systematic_resample <- function (weights)
  .Call(P_systematic_resampling,weights)
