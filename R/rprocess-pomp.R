## simulate the process model
setMethod(
          'rprocess',
          'pomp',
          function (object, xstart, times, params, ...) { # the package algorithms will only use these arguments
            ntimes <- length(times)
            if (ntimes<2)
              stop("rprocess error: no transitions: no work to do",call.=FALSE)
            if (ncol(params)!=ncol(xstart))
              stop("rprocess error: number of columns of 'params' and 'xstart' do not agree",call.=FALSE)
            x <- try(
                     do.call(
                             object@rprocess,
                             c(
                               list(
                                    xstart=xstart,
                                    times=times,
                                    params=params
                                    ),
                               object@userdata   # the userdata gets sent as extra arguments to the user's rprocess function
                               )
                             ),
                     silent=FALSE
                     )
            if (inherits(x,'try-error'))
              stop("rprocess error: error in user 'rprocess'",call.=FALSE)
            dim.x <- dim(x)
            if (length(dim.x)!=3 || any(dim.x!=c(dim(xstart),ntimes)) || is.null(rownames(x)))
              stop(
                   "rprocess error: user 'rprocess' must return an array of dimensions ",
                   nrow(xstart)," x ",ncol(xstart)," x ",ntimes," with rownames",
                   call.=FALSE
                   )
            x
          }
          )
