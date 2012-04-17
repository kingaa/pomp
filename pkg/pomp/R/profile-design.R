profileDesign <- function (..., lower, upper, nprof) {
  prof <- list(...)
  pvars <- names(prof)
  if (any(pvars==""))
    stop(sQuote("profileDesign"),": you cannot profile over an unnamed variable!")
  ovars <- names(lower)
  if (!all(sort(ovars)==sort(names(upper))))
    stop(sQuote("profileDesign"),": names of ",sQuote("lower")," and ",sQuote("upper")," must match!")
  vars <- ovars[!(ovars%in%pvars)]
  x <- as.matrix(expand.grid(...))
  y <- as.matrix(sobolDesign(lower=lower,upper=upper,nseq=nprof))
  z <- array(dim=c(nrow(x)*nrow(y),ncol(x)+ncol(y)))
  colnames(z) <- c(colnames(x),colnames(y))
  i <- 1
  for (j in seq_len(nrow(x))) {
    for (k in seq_len(nrow(y))) {
      z[i,] <- c(x[j,],y[k,])
      i <- i+1
    }
  }
  as.data.frame(z)
}
