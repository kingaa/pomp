## evaluate the measurement model density function

dmeasure.internal <- function (object, y, x, times, params, log = FALSE, .getnativesymbolinfo = TRUE, ...) {
  .Call(do_dmeasure,object,y,x,times,params,log,.getnativesymbolinfo)
}

setMethod("dmeasure","pomp",
          function (object, y, x, times, params, log = FALSE, ...)
          dmeasure.internal(object=object,y=y,x=x,times=times,params=params,log=log,...)
          )
