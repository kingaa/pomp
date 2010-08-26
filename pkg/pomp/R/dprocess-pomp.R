dprocess <- function (object, x, times, params, log = FALSE, ...)
  stop("function ",sQuote("dprocess")," is undefined for objects of class ",sQuote(class(object)))
setGeneric('dprocess')  

## evaluate the process model density function
setMethod(
          'dprocess',
          'pomp',
## SEXP do_dprocess (SEXP object, SEXP x, SEXP times, SEXP params, SEXP log)
          function (object, x, times, params, log = FALSE, ...) {
            x <- try(
                     .Call(do_dprocess,object,x,times,params,log),
                     silent=FALSE
                     )
            if (inherits(x,'try-error'))
              stop("dprocess error: error in user ",sQuote("dprocess"),call.=FALSE)
            x
          }
          )
