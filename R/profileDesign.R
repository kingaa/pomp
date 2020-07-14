##' @description
##' Obsolete: replaced by \code{\link{profile_design}}. 
##' 
##' @rdname deprecated
##' @return 
##' \code{profileDesign} returns a data frame with \code{nprof} points per profile point.
##' @param nprof The number of points per profile point.
##' @param type the type of design to use.
##' \code{type="sobol"} uses \code{\link{sobolDesign}};
##' \code{type="runif"} uses \code{\link{runifDesign}}.
##' @param stringsAsFactors should character vectors be converted to factors?
##' 
##' @export
profileDesign <- function (...,
  lower, upper, nprof,
  type = c("sobol","runif"),
  stringsAsFactors = getOption("stringsAsFactors",FALSE)
) {
  .Deprecated("profile_design",package="pomp")
  ep <- "profileDesign"
  type <- match.arg(type)
  prof <- list(...)
  pvars <- names(prof)
  if (any(pvars==""))
    pStop(ep,"you cannot profile over an unnamed variable!")
  ovars <- names(lower)
  if (!all(sort(ovars)==sort(names(upper))))
    pStop(ep,"names of ",sQuote("lower")," and ",sQuote("upper")," must match!")
  x <- expand.grid(...,stringsAsFactors=stringsAsFactors)
  z <- vector(mode='list',length=nrow(x))
  y <- switch(type,
    runif=runifDesign(lower=lower,upper=upper,nseq=nprof),
    sobolDesign(lower=lower,upper=upper,nseq=nprof)
  )
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
