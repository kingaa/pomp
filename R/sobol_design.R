##' @description
##' \code{sobol_design} generates a Latin hypercube design based on the Sobol' low-discrepancy sequence.
##' @rdname design
##' 
##' @details
##' The Sobol' sequence generation is performed using codes from the
##' \href{http://ab-initio.mit.edu/nlopt/}{\pkg{NLopt} library} by S. Johnson.
##' @return
##' \code{sobol_design} returns a data frame with \code{nseq} rows and one column for each variable named in \code{lower} and \code{upper}.
##' @param lower,upper named numeric vectors giving the lower and upper bounds
##' of the ranges, respectively.
##' @param nseq Total number of points requested.
##' 
##' @export
sobol_design <- function (lower = numeric(0), upper = numeric(0), nseq) {
  ep <- "sobol_design"
  if (length(lower)!=length(upper))
    pStop(ep,sQuote("lower")," and ",sQuote("upper")," must have same length.")
  lnames <- names(lower)
  if (is.null(lnames))
    pStop(ep,sQuote("lower")," and ",sQuote("upper")," must be named vectors.")
  if (!all(sort(lnames)==sort(names(upper))))
    pStop(ep,"names of ",sQuote("lower")," and ",sQuote("upper")," must match.")
  upper <- upper[lnames]
  ranges <- lapply(seq_along(lnames),function(k)c(lower[k],upper[k]))
  names(ranges) <- lnames
  tryCatch(
    sobol(ranges,n=as.integer(nseq)),
    error = function (e) pStop(ep,conditionMessage(e))
  )
}

sobol <- function (vars, n) {
  d <- length(vars)
  if (!is.finite(n) || (n > 1073741824L))
    pStop_("too many points requested.")
  x <- .Call(P_sobol_sequence,as.integer(d),as.integer(n))
  y <- vapply(
    seq_len(d),
    function (k) {
      vars[[k]][1L]+(vars[[k]][2L]-vars[[k]][1L])*x[k,]
    },
    numeric(n)
  )
  colnames(y) <- names(vars)
  as.data.frame(y)
}
