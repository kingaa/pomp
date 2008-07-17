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
              stop("rmeasure error: error in user 'rmeasure'",call.=FALSE)
            val
          }
          )
