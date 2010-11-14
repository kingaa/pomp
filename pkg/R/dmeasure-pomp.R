setGeneric("dmeasure",function(object,...)standardGeneric("dmeasure"))

## evaluate the measurement model density function
setMethod(
          'dmeasure',
          'pomp',
          function (object, y, x, times, params, log = FALSE, ...)
            .Call(do_dmeasure,object,y,x,times,params,log),
          )
