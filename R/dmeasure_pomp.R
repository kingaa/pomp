## evaluate the measurement model density function

setMethod(
  "dmeasure",
  signature=signature(object="pomp"),
  definition=function (object, y, x, times, params, log = FALSE, ...) {
    dmeasure.internal(object=object,y=y,x=x,times=times,
      params=params,log=log,...)
  }
)

dmeasure.internal <- function (object, y, x, times, params, log = FALSE,
  .getnativesymbolinfo = TRUE, ...) {
  tryCatch(
    {
      storage.mode(y) <- "double"
      storage.mode(x) <- "double"
      storage.mode(params) <- "double"
      pompLoad(object)
      rv <- .Call(do_dmeasure,object,y,x,times,params,log,.getnativesymbolinfo)
      pompUnload(object)
      rv
    },
    error = function (e) {
      stop("in ",sQuote("dmeasure"),": ",conditionMessage(e),call.=FALSE)
    }
  )
}
