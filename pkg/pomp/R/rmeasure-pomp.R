## simulate the measurement model

setGeneric("rmeasure",function(object,...)standardGeneric("rmeasure"))

rmeasure.internal <- function (object, x, times, params, .getnativesymbolinfo = TRUE, ...) {
  .Call(do_rmeasure,object,x,times,params,.getnativesymbolinfo)
}

setMethod("rmeasure","pomp",
          function (object, x, times, params, ...)
          rmeasure.internal(object=object,x=x,times=times,params=params,...)
          )
