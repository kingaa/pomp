##' @description
##' \code{slice_design} generates points along slices through a specified point.
##' @rdname design
##' @return
##' \code{slice_design} returns a data frame with one row per point.
##' The \sQuote{slice} variable indicates which slice the point belongs to.
##' @param center \code{center} is a named numeric vector specifying the point
##' through which the slice(s) is (are) to be taken.
##'
##' @export
slice_design <- function (center, ...) {
  slices <- list(...)
  if ((!is.numeric(center))||is.null(names(center)))
    pStop(sQuote("center")," must be a named numeric vector")
  slnm <- names(slices)
  if (any(slnm==""))
    pStop("cannot slice along an unnamed parameter.")
  if (!all(slnm%in%names(center))) {
    problems <- slnm[!(slnm%in%names(center))]
    pStop(
      ngettext(length(problems),"variable ","variables "),
      paste(lapply(problems,sQuote),collapse=","),
      ngettext(length(problems)," does "," do "),
      "not appear in ",sQuote("center")
    )
  }
  nslice <- length(slices)
  nvars <- length(center)

  y <- vector(mode="list",length=nslice)
  for (k in seq_len(nslice)) {
    y[[k]] <- as.data.frame(
      matrix(
        data=center,
        nrow=length(slices[[k]]),
        ncol=nvars,
        byrow=TRUE,
        dimnames=list(names(slices[[k]]),names(center))
      )
    )
    y[[k]][[slnm[k]]] <- as.numeric(slices[[k]])
    y[[k]]$slice <- slnm[k]
  }
  y <- do.call(rbind,y)
  y$slice <- as.factor(y$slice)
  y
}
