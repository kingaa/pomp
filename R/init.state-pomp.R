## initialize the process model
setMethod(
          'init.state',
          'pomp',
          function (object, params, t0, ...) { # the package algorithms will only use these arguments
            if (missing(t0))
              t0 <- object@t0
            nreps <- NCOL(params)
            x <- try(
                     apply(
                           as.matrix(params),
                           2,
                           function(p).Call(do_init_state,object,p,t0)
                           ),
                     silent=FALSE
                     )
            if (inherits(x,'try-error'))
              stop("init.state error: error in user 'initializer'",call.=FALSE)
            x
          }
          )
