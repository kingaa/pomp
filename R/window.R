##' Window
##'
##' Restrict to a portion of a time series.
##'
##' @name window
##' @docType methods
##' @rdname window
##' @include pomp_class.R
NULL

##' @importFrom stats window
setGeneric("window")

##' @rdname window
##' @param x a \sQuote{pomp} object or object of class extending \sQuote{pomp}
##' @param start,end the left and right ends of the window, in units of time
##' @param \dots ignored
##' @export
setMethod(
  "window",
  signature=signature(x="pomp"),
  definition=function (x, start, end, ...) {
    tm <- time(x,t0=FALSE)
    if (missing(start)) start <- tm[1L]
    if (missing(end)) end <- tm[length(tm)]
    if (!(is.numeric(start) && is.finite(start) && length(start)==1 &&
        is.numeric(end) && is.finite(end) && length(end)==1))
      pStop(who="window",sQuote("start")," and ",sQuote("end"),
        " must be finite times.")
    if (!isTRUE(start <= end))
      pStop(who="window",sQuote("start"),
        " must not be later than ",sQuote("end"),".")
    tm <- tm[(tm>=start)&(tm<=end)]
    time(x,t0=FALSE) <- tm
    x
  }
)
