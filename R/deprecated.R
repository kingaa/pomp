conv.rec <- function (object, ...) {
  warning(sQuote("conv.rec")," is deprecated and will be removed in a ",
    "forthcoming release.  Please use ",sQuote("traces")," instead.",
    call.=FALSE)
  traces(object,...)
}

values <- function (object, ...) {
  warning(sQuote("values")," is deprecated and will be removed in a ",
    "forthcoming release.  Please use ",sQuote("probevals")," or ",
    sQuote("as.data.frame")," instead.",call.=FALSE)
  as(object,"data.frame")
}
