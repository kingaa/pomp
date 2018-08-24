##' @rdname design
##' @export
sobolDesign <- function (lower = numeric(0), upper = numeric(0), nseq) {
  ep <- "sobolDesign"
  if (length(lower)!=length(upper))
    pStop(ep,sQuote("lower")," and ",sQuote("upper")," must have same length.")
  lnames <- names(lower)
  if (is.null(lnames))
    pStop(ep,sQuote("lower")," and ",sQuote("upper")," must be named vectors.")
  if (!all(sort(lnames)==sort(names(upper))))
    pStop(ep,"names of ",sQuote("lower")," and ",sQuote("upper")," must match.")
  upper <- upper[lnames]
  ranges <- lapply(lnames,function(k)c(lower[k],upper[k]))
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
  x <- .Call(sobol_sequence,as.integer(d),as.integer(n))
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
