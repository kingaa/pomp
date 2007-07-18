## simulate the measurement model
setMethod(
          'rmeasure',
          'pomp',
          function (object, x, times, params, ...) {
            if (length(dim(x))!=3)
              stop("rmeasure error: 'x' must be a rank-3 array")
            ntimes <- length(times)
            nsims <- dim(x)[2]
            if (ntimes!=dim(x)[3])
              stop("rmeasure error: length of 'times' does not agree with dim(x)[3]")
            if (nsims!=ncol(params))
              stop("rmeasure error: number of columns of 'params' and 'x' do not agree")
            y <- try(
                     do.call(
                             object@rmeasure,
                             c(
                               list(
                                    x=x,
                                    times=times,
                                    params=params
                                    ),
                               object@userdata
                               )
                             ),
                     silent=T
                     )
            if (inherits(y,'try-error'))
              stop("rmeasure error: error in user 'rmeasure'\n",y)
            nobs <- nrow(object@data)
            dim.y <- dim(y)
            if (length(dim.y)!=3 || any(dim.y!=c(nobs,nsims,ntimes)))
              stop("rmeasure error: user 'rmeasure' must return an array of dimensions ",
                   nobs,"x",nsims,"x",ntimes)
            rownames(y) <- rownames(object@data)
            y
          }
          )
