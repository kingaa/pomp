##' @rdname design
##' @export
sliceDesign <- function (center, ...) {
  slices <- list(...)
  if ((!is.numeric(center))||is.null(names(center)))
    pStop("sliceDesign",sQuote("center")," must be a named numeric vector")
  slnm <- names(slices)
  if (any(slnm==""))
    pStop("sliceDesign","cannot slice along an unnamed parameter.")
  if (!all(slnm%in%names(center))) {
    problems <- slnm[!(slnm%in%names(center))]
    pStop("sliceDesign",
      ngettext(length(problems),"variable ","variables "),
      paste(sapply(problems,sQuote),collapse=","),
      ngettext(length(problems)," does "," do "),
      "not appear in ",sQuote("center"))
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
