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
              stop("skeleton error: error in user 'skeleton'",call.=FALSE)
            x
          }
          )
