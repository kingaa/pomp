##' @description
##' \code{runif_design} generates a design based on random samples from a multivariate uniform distribution.
##' @rdname design
##' @return
##' \code{runif_design} returns a data frame with \code{nseq} rows and one column for each variable named in \code{lower} and \code{upper}.
##'
##' @export
runif_design <- function (lower = numeric(0), upper = numeric(0), nseq) {
  if (length(lower)!=length(upper))
    pStop(sQuote("lower")," and ",sQuote("upper")," must have same length.")
  lnames <- names(lower)
  if (is.null(lnames))
    pStop(sQuote("lower")," and ",sQuote("upper")," must be named vectors.")
  if (!all(sort(lnames)==sort(names(upper))))
    pStop("names of ",sQuote("lower")," and ",sQuote("upper")," must match.")
  upper <- upper[lnames]
  if (!all(upper>=lower))
    pStop("upper values should be at least as large as lower ones.")
  nseq <- as.integer(nseq)
  if (nseq < 0)
    pStop(sQuote("nseq"),"< 0.")
  y <- matrix(
    data=runif(n=nseq*length(lower),min=lower,max=upper),
    nrow=nseq,ncol=length(lower),
    byrow=TRUE
  )
  colnames(y) <- lnames
  as.data.frame(y)
}
