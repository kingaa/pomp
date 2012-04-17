setGeneric("dprocess",function(object,...)standardGeneric("dprocess"))

## evaluate the process model density function
setMethod(
          'dprocess',
          'pomp',
          function (object, x, times, params, log = FALSE, ...)
            .Call(do_dprocess,object,x,times,params,log)
          )
