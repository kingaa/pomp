## simulate the process model
setMethod(
          'rprocess',
          'pomp',
          function (object, xstart, times, params, ...) { # the package algorithms will only use these arguments
            x <- try(
                     .Call(do_rprocess,object,xstart,times,params),
                     silent=FALSE
                     )
            if (inherits(x,'try-error'))
              stop("rprocess error: error in user 'rprocess'",call.=FALSE)
            x
          }
          )
