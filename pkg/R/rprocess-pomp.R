setGeneric("rprocess",function(object,...)standardGeneric("rprocess"))

## simulate the process model
setMethod(
          'rprocess',
          'pomp',
          function (object, xstart, times, params, offset = 0, ...)
            .Call(do_rprocess,object,xstart,times,params,offset)
          )
