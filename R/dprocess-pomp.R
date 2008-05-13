## evaluate the process model density function
setMethod(
          'dprocess',
          'pomp',
          function (object, x, times, params, log = FALSE, ...) {
            if (length(dim(x))!=3)
              stop("dprocess error: 'x' must be a rank-3 array",call.=FALSE)
            if (length(dim(params))!=2)
              stop("dprocess error: 'params' must be a rank-2 array",call.=FALSE)
            ntimes <- length(times)
            if (length(times)!=dim(x)[3])
              stop("dprocess error: dim(x)[3] != length(times)",call.=FALSE)
            nsims <- dim(x)[2]
            if (nsims!=NCOL(params))
              stop("dprocess error: number of columns of 'params' and 'x' do not agree",call.=FALSE)
            d <- try(
                     do.call(
                             object@dprocess,
                             c(
                               list(
                                    x=x,
                                    times=times,
                                    params=params,
                                    log=log
                                    ),
                               object@userdata
                               )
                             ),
                     silent=FALSE
                     )
            if (inherits(d,'try-error'))
              stop("dprocess error: error in user 'dprocess'",call.=FALSE)
            if (any(dim(d)!=c(nsims,ntimes-1)))
              stop(
                   "dprocess error: user 'dprocess' must return a matrix with one value per replicate per transition",
                   call.=FALSE
                   )
            d
          }
          )

