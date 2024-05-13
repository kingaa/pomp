##' Lookup table
##'
##' Interpolate values from a lookup table
##'
##' @rdname lookup
##' @name lookup
##' @family interpolation
##' @include covariate_table.R
##' @param table a \sQuote{covartable} object created by a call to \code{\link{covariate_table}}
##' @param t numeric vector; times at which interpolated values of the covariates in \code{table} are required.
##' @return
##' A numeric vector or matrix of the interpolated values.
##' @inheritSection covariates Extrapolation
##' @export

lookup <- function (table, t) {
  d <- .Call(P_lookup_in_table,table,t)
  data.frame(t=t,t(d))
}

##' @rdname covariate_table
##' @return
##' \code{repair_lookup_table} returns a lookup table with entries at the provided values of \code{t}.
##' @param table a \sQuote{covartable} object created by a call to \code{\link{covariate_table}}
##' @param t numeric vector;
##' times at which interpolated values of the covariates in \code{table} are required.
##' @details
##' \code{repair_lookup_table} applies \code{\link{lookup}} at the provided values of \code{t} and returns the resulting lookup table.
##' \strong{\code{repair_lookup_table} should be considered experimental: its interface may change without notice}.
##' @export
repair_lookup_table <- function (table, t) {
  covariate_table(
    lookup(table,t=t),
    order=if (table@order==0L) "constant" else "linear",
    times="t"
  )
}
