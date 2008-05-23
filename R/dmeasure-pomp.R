## evaluate the measurement model density function
setMethod(
          'dmeasure',
          'pomp',
          function (object, y, x, times, params, log = FALSE, ...) {
            x <- try(
                     .Call(do_dmeasure,object,y,x,times,params,log),
                     silent=TRUE
                     )
            if (inherits(x,'try-error'))
              stop("dmeasure error: error in user 'dmeasure'",call.=FALSE)
            x
          }
          )
