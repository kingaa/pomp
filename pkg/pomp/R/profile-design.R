profileDesign <- function (..., lower, upper, nprof,
                           stringsAsFactors = default.stringsAsFactors()) {
  prof <- list(...)
  pvars <- names(prof)
  if (any(pvars==""))
    stop(sQuote("profileDesign"),": you cannot profile over an unnamed variable!")
  ovars <- names(lower)
  if (!all(sort(ovars)==sort(names(upper))))
    stop(sQuote("profileDesign"),": names of ",sQuote("lower")," and ",sQuote("upper")," must match!")
  x <- expand.grid(...,stringsAsFactors=stringsAsFactors)
  y <- sobolDesign(lower=lower,upper=upper,nseq=nprof)
  z <- vector(mode='list',length=nrow(x)*nprof)
  for (i in seq_len(nrow(x))) {
    z[[i]] <- data.frame(
                         x[i,,drop=FALSE],y,
                         check.rows=FALSE,
                         check.names=FALSE,
                         row.names=seq_len(nprof)
                         )
  }
  do.call(rbind,z)
}
