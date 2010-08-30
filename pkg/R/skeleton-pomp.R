skeleton <- function (object, x, t, params, ...)
  stop("function ",sQuote("skeleton")," is undefined for objects of class ",sQuote(class(object)))
setGeneric('skeleton')  

## evaluate the measurement model density function
setMethod(
          'skeleton',
          'pomp',
          function (object, x, t, params, ...) {
            x <- try(
                     .Call(do_skeleton,object,x,t,params),
                     silent=FALSE
                     )
            if (inherits(x,'try-error'))
              stop("skeleton error: error in user ",sQuote("skeleton"),call.=FALSE)
            x
          }
          )
