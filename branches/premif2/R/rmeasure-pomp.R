## simulate the measurement model

setGeneric("rmeasure",function(object,...)standardGeneric("rmeasure"))

rmeasure.internal <- function (object, x, times, params, ...) {
  fun <- get.pomp.fun(object@rmeasure)
  .Call(do_rmeasure,object,x,times,params,fun)
}

setMethod("rmeasure","pomp",rmeasure.internal)
