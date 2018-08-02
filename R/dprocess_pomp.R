## evaluate the process model density function

dprocess.internal <- function (object, x, times, params, log = FALSE, .getnativesymbolinfo = TRUE, ...) {
    storage.mode(x) <- "double"
    storage.mode(params) <- "double"
    pompLoad(object)
    on.exit(pompUnload(object))
    rv <- tryCatch(
      .Call(do_dprocess,object,x,times,params,log,.getnativesymbolinfo),
      error=function (e) {
        stop("in ",sQuote("dprocess"),": ",conditionMessage(e),call.=FALSE)
      })
    rv
}

setMethod(
    "dprocess",
    signature=signature(object="pomp"),
    definition = function (object, x, times, params, log = FALSE, ...)
        dprocess.internal(object=object,x=x,times=times,params=params,log=log,...)
)
