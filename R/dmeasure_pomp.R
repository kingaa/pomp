## evaluate the measurement model density function

dmeasure.internal <- function (object, y, x, times, params, log = FALSE, .getnativesymbolinfo = TRUE, ...) {
  pompLoad(object)
  rv <- .Call(do_dmeasure,object,y,x,times,params,log,.getnativesymbolinfo)
  pompUnload(object)
  rv
}

setMethod(
    "dmeasure",
    signature=signature(object="pomp"),
    definition=function (object, y, x, times, params, log = FALSE, ...)
        dmeasure.internal(object=object,y=y,x=x,times=times,params=params,log=log,...)
)
