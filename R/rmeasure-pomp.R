## simulate the measurement model
setMethod(
          'rmeasure',
          'pomp',
          function (object, x, times, params, ...) {
            x <- try(
                     .Call(do_rmeasure,object,x,times,params),
                     silent=FALSE
                     )
            if (inherits(x,'try-error'))
              stop("rmeasure error: error in user 'rmeasure'",call.=FALSE)
            x
          }
          )
