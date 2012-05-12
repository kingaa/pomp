## simulate the process model

setGeneric("rprocess",function(object,...)standardGeneric("rprocess"))

rprocess.internal <- function (object, xstart, times, params, offset = 0, ...)
  .Call(do_rprocess,object,xstart,times,params,offset)

setMethod("rprocess","pomp",rprocess.internal)
