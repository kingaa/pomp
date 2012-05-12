## evaluate the measurement model density function
setGeneric("dmeasure",function(object,...)standardGeneric("dmeasure"))

dmeasure.internal <- function (object, y, x, times, params, log = FALSE, ...) {
  fun <- get.pomp.fun(object@dmeasure)
  .Call(do_dmeasure,object,y,x,times,params,log,fun)
}

setMethod("dmeasure","pomp",dmeasure.internal)
