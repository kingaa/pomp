## evaluate the process model density function
setGeneric("dprocess",function(object,...)standardGeneric("dprocess"))

dprocess.internal <- function (object, x, times, params, log = FALSE, ...)
  .Call(do_dprocess,object,x,times,params,log)

setMethod("dprocess","pomp",dprocess.internal)
