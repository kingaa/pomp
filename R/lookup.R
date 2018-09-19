##' Lookup table
##'
##' Interpolate values from a lookup table
##'
##' @rdname lookup
##' @name lookup
##' @keywords internal
##' @concept interpolation
##' @seealso covariate_table
##'
##' @param table a \sQuote{covartable} object created by a call to \code{\link{covariate_table}}
##' @param t numeric vector; times at which interpolated values of the covariates in \code{table} are required.
##'
##' @return
##' A numeric vector or matrix of the interpolated values.
##'
##' @details
##' A warning will be generated if extrapolation is performed.
##' @export
##'
lookup <- function (table, t) {
  d <- .Call(P_lookup_in_table,table,t)
  data.frame(t=t,t(d))
}
