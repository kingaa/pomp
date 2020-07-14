##' @description
##' Obsolete: replaced by \code{sobol_design}.
##' @rdname deprecated
##' 
##' @details
##' The Sobol' sequence generation is performed using codes from the
##' \href{http://ab-initio.mit.edu/nlopt/}{\pkg{NLopt} library} by S. Johnson.
##' @return
##' \code{sobolDesign} returns a data frame with \code{nseq} rows and one column for each variable named in \code{lower} and \code{upper}.
##' @param lower,upper named numeric vectors giving the lower and upper bounds
##' of the ranges, respectively.
##' @param nseq Total number of points requested.
##' 
##' @export
sobolDesign <- function (lower = numeric(0), upper = numeric(0), nseq) {
  .Deprecated("sobol_design",package="pomp")
  ep <- "sobolDesign"
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
