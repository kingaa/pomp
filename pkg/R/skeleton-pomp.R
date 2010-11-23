skeleton <- function (object, x, t, params, ...)
  stop("function ",sQuote("skeleton")," is undefined for objects of class ",sQuote(class(object)))

setGeneric('skeleton')

## evaluate the measurement model density function
setMethod(
          'skeleton',
          'pomp',
          function (object, x, t, params, ...) {
            skel <- .Call(get_pomp_fun,object@skeleton)
            .Call(do_skeleton,object,x,t,params,skel)
          }
          )
