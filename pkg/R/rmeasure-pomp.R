setGeneric("rmeasure",function(object,...)standardGeneric("rmeasure"))

## simulate the measurement model
setMethod(
          'rmeasure',
          'pomp',
          function (object, x, times, params, ...)
            .Call(do_rmeasure,object,x,times,params),
          )
