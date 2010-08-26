rmeasure <- function (object, x, times, params, ...)
  stop("function ",sQuote("rmeasure")," is undefined for objects of class ",sQuote(class(object)))
setGeneric('rmeasure')  

## simulate the measurement model
setMethod(
          'rmeasure',
          'pomp',
          function (object, x, times, params, ...) {
            val <- try(
                       .Call(do_rmeasure,object,x,times,params),
                       silent=FALSE
                       )
            if (inherits(val,'try-error'))
              stop("rmeasure error: error in user ",sQuote("rmeasure"),call.=FALSE)
            val
          }
          )
