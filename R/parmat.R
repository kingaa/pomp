##' Create a matrix of parameters
##'
##' \code{parmat} is a utility that makes a vector of parameters suitable for
##' use in \pkg{pomp} functions.
##'
##'
##' @param params named numeric vector or matrix of parameters.
##' @param nrep number of replicates (columns) desired.
##' @param names optional character; column names.
##' @return \code{parmat} returns a matrix consisting of \code{nrep} copies of
##' \code{params}.
##' @author Aaron A. King
##' @example examples/ricker-bifdiag.R
##' @export
parmat <- function (params, nrep = 1, names = NULL) {
  d <- dim(params)
  if (is.null(d) || length(d) == 1) {
    matrix(data=params,nrow=length(params),ncol=nrep,
      dimnames=list(variable=names(params),rep=names))
  } else if (length(d) == 2) {
    matrix(data=params,nrow=nrow(params),ncol=ncol(params)*nrep,
      dimnames=list(variable=rownames(params),rep=names))
  } else {
    matrix(data=params,nrow=nrow(params),ncol=prod(d[-1])*nrep,
      dimnames=list(variable=rownames(params),rep=names))
  }
}
