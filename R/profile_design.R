##' @description
##' \code{profile_design} generates a data-frame where each row can be used as the starting point for a profile likelihood calculation.
##' @concept profile likelihood
##' 
##' @rdname design
##' @return 
##' \code{profile_design} returns a data frame with \code{nprof} points per profile point.
##' @param nprof The number of points per profile point.
##' @param type the type of design to use.
##' \code{type="runif"} uses \code{\link{runif_design}}.
##' \code{type="sobol"} uses \code{\link{sobol_design}};
##' @param stringsAsFactors should character vectors be converted to factors?
##' 
##' @export
profile_design <- function (...,
  lower, upper, nprof,
  type = c("runif","sobol"),
  stringsAsFactors = getOption("stringsAsFactors",FALSE)
) {
  ep <- "profile_design"
  type <- match.arg(type)
  prof <- list(...)
  pvars <- names(prof)
  if (any(pvars==""))
    pStop(ep,"you cannot profile over an unnamed variable!")
  ovars <- names(lower)
  if (!all(sort(ovars)==sort(names(upper))))
    pStop(ep,"names of ",sQuote("lower")," and ",sQuote("upper")," must match!")
  x <- expand.grid(...,stringsAsFactors=stringsAsFactors)
  n <- nrow(x)
  z <- vector(mode='list',length=n)
  y <- switch(type,
    runif=runif_design(lower=lower,upper=upper,nseq=n*nprof),
    sobol_design(lower=lower,upper=upper,nseq=n*nprof)
  )
  y <- split(y,rep(seq_len(n),each=nprof))
  for (i in seq_len(n)) {
    z[[i]] <- data.frame(
      x[i,,drop=FALSE], y[[i]],
      check.rows=FALSE,
      check.names=FALSE,
      row.names=seq_len(nprof)
    )
  }
  do.call(rbind,z)
}
