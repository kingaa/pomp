sobol <- function (vars, n) {
  if (!is.list(vars) || is.null(names(vars)))
    stop("sobol error: ",sQuote("vars")," must be a named list")
  if (!all(sapply(vars,function(x)is.numeric(x)&&(length(x)==2))))
    stop("sobol error: each entry in ",sQuote("vars")," must specify a range")
  d <- length(vars)
  x <- .Call("sobol_sequence",as.integer(c(d,n)))
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

sobolDesign <- function (lower = numeric(0), upper = numeric(0), nseq) {
  if (length(lower)!=length(upper))
    stop(sQuote("lower")," and ",sQuote("upper")," must have same length")
  lnames <- names(lower)
  if (is.null(lnames))
    stop(sQuote("lower")," and ",sQuote("upper")," must be named vectors")
  if (!all(sort(lnames)==sort(names(upper))))
    stop("names of ",sQuote("lower")," and ",sQuote("upper")," must match")
  upper <- upper[lnames]
  ranges <- lapply(lnames,function(k)c(lower[k],upper[k]))
  names(ranges) <- lnames
  sobol(ranges,n=as.integer(nseq))
}
