## simulate the measurement model

setMethod(
  "rmeasure",
  signature=signature(object="pomp"),
  definition=function (object, x, times, params, ...)
    rmeasure.internal(object=object,x=x,times=times,params=params,...)
)

rmeasure.internal <- function (object, x, times, params,
  .getnativesymbolinfo = TRUE, ...) {
  tryCatch(
    {
      storage.mode(x) <- "double"
      storage.mode(params) <- "double"
      pompLoad(object)
      rv <- .Call(do_rmeasure,object,x,times,params,.getnativesymbolinfo)
      pompUnload(object)
      rv
    },
    error = function (e) {
      stop("in ",sQuote("rmeasure"),": ",conditionMessage(e),call.=FALSE)
    }
  )
}
