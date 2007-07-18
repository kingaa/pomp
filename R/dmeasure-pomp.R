## evaluate the measurement model density function
setMethod(
          'dmeasure',
          'pomp',
          function (object, y, x, times, params, log = FALSE, ...) {
            if (length(dim(y))!=2)
              stop("dmeasure error: 'y' must be a rank-2 array")
            if (length(dim(x))!=3)
              stop("dmeasure error: 'x' must be a rank-3 array")
            if (dim(y)[2]!=dim(x)[3])
              stop("dmeasure error: dim(x)[3] != dim(y)[2]")
            ntimes <- length(times)
            if (length(times)!=dim(y)[2])
              stop("dmeasure error: dim(y)[2] != length(times)")
            nsims <- ncol(x)
            if (nsims!=NCOL(params))
              stop("dmeasure error: number of columns of 'params' and 'x' do not agree")
            d <- try(
                     do.call(
                             object@dmeasure,
                             c(
                               list(
                                    y=y,
                                    x=x,
                                    times=times,
                                    params=params,
                                    log=log
                                    ),
                               object@userdata
                               )
                             ),
                     silent=T
                     )
            if (inherits(d,'try-error'))
              stop("dmeasure error: error in user 'dmeasure'\n",d)
            if (any(dim(d)!=c(nsims,ntimes)))
              stop("dmeasure error: user 'dmeasure' must return one value per replicate per time")
            d
          }
          )

