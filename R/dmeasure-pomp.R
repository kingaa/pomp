## evaluate the measurement model density function
setMethod(
          'dmeasure',
          'pomp',
          function (object, y, x, times, params, log = FALSE, ...) {
            val <- try(
                       .Call(do_dmeasure,object,y,x,times,params,log),
                       silent=FALSE
                       )
            if (inherits(val,'try-error'))
              stop("dmeasure error: error in user 'dmeasure'",call.=FALSE)
            val
          }
          )
