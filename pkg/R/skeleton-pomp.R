setGeneric("skeleton",function(object,...)standardGeneric("skeleton"))

## evaluate the measurement model density function
setMethod(
          'skeleton',
          'pomp',
          function (object, x, t, params, ...)
            .Call(do_skeleton,object,x,t,params),
          )
