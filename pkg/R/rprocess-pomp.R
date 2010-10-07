rprocess <- function (object, xstart, times, params, ...)
  stop("function ",sQuote("rprocess")," is undefined for objects of class ",sQuote(class(object)))
setGeneric('rprocess')  

## simulate the process model
setMethod(
          'rprocess',
          'pomp',
          function (object, xstart, times, params, offset = 0, ...)
            .Call(do_rprocess,object,xstart,times,params,offset)
          )
