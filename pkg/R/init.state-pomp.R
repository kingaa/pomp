## initialize the process model
setMethod(
          'init.state',
          'pomp',
          function (object, params, t0, ...) { # the package algorithms will only use these arguments
            if (missing(t0))
              t0 <- object@t0
            x <- try(
                     .Call(do_init_state,object,as.matrix(params),t0),
                     silent=FALSE
                     )
            if (inherits(x,'try-error'))
              stop("init.state error: error in user ",sQuote("initializer"),call.=FALSE)
            x
          }
          )
