## evaluate the process model density function
setGeneric("dprocess",function(object,...)standardGeneric("dprocess"))

dprocess.internal <- function (object, x, times, params, log = FALSE, .getnativesymbolinfo = TRUE, ...)
  .Call(do_dprocess,object,x,times,params,log,.getnativesymbolinfo)

setMethod("dprocess","pomp",
          function (object, x, times, params, log = FALSE, ...)
          dprocess.internal(object=object,x=x,times=times,params=params,log=log,...)
          )
