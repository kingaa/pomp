conv.rec <- function (object, ...) {
  warning(sQuote("conv.rec")," is deprecated and will be removed in a ",
    "forthcoming release. Please use ",sQuote("traces")," instead.",
    call.=FALSE)
  traces(object,...)
}
