## initialize the process model
setMethod(
          'init.state',
          'pomp',
          function (object, params, t0, ...) { # the package algorithms will only use these arguments
            if (missing(t0))
              t0 <- object@t0
            x <- vector(mode='list',length=NCOL(params))
            params <- as.matrix(params)
            for (k in seq(NCOL(params))) {
              x[[k]] <- try(
                            do.call(
                                    object@initializer,
                                    c(
                                      list(
                                           params=params[,k],
                                           t0=t0
                                           ),
                                      object@userdata   # the userdata gets sent as extra arguments to the user's rprocess function
                                      )
                                    ),
                            silent=FALSE
                            )
              if (inherits(x[[k]],'try-error'))
                stop("init.state error: error in user 'initializer'")
              if ((!is.numeric(x[[k]]))||(is.array(x[[k]]))||is.null(names(x[[k]])))
                stop("'initializer' must return a named vector")
            }
            do.call(cbind,x)
          }
          )
