sliceDesign <- function (center, ...) {
  slices <- list(...)
  if ((!is.numeric(center))||is.null(names(center)))
    stop(sQuote("sliceDesign"),": ",sQuote("center")," must be a named numeric vector")
  slnm <- names(slices)
  if (any(slnm==""))
    stop(sQuote("sliceDesign"),": you cannot slice along an unnamed parameter")
  if (!all(slnm%in%names(center))) {
    problems <- slnm[!(slnm%in%names(center))]
    stop(
         sQuote("sliceDesign"),
         " error: variable(s) ",
         sQuote(paste(problems,collapse=",")),
         " do not appear in ",
         sQuote("center")
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
