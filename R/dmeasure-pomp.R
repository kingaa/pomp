## evaluate the measurement model density function
setMethod(
          'dmeasure',
          'pomp',
          function (object, y, x, times, params, log = FALSE, ...) {
            if (length(dim(y))!=2)
              stop("dmeasure error: 'y' must be a rank-2 array",call.=FALSE)
            if (length(dim(x))!=3)
              stop("dmeasure error: 'x' must be a rank-3 array",call.=FALSE)
            if (length(dim(params))!=2)
              stop("dmeasure error: 'params' must be a rank-2 array",call.=FALSE)
            if (dim(params)[2]!=dim(x)[2])
              stop("dmeasure error: ncol(x) != ncol(params)",call.=FALSE)
            ntimes <- length(times)
            if (length(times)!=dim(y)[2])
              stop("dmeasure error: ncol(y) != length(times)",call.=FALSE)
            if (length(times)!=dim(x)[3])
              stop("dmeasure error: dim(x)[3] != length(times)",call.=FALSE)
            nsims <- ncol(x)
            if (nsims!=NCOL(params))
              stop("dmeasure error: number of columns of 'params' and 'x' do not agree",call.=FALSE)
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
                     silent=FALSE
                     )
            if (inherits(d,'try-error'))
              stop("dmeasure error: error in user 'dmeasure'",call.=FALSE)
            if (any(dim(d)!=c(nsims,ntimes)))
              stop(
                   "dmeasure error: user 'dmeasure' must return a matrix with one value per replicate per time",
                   call.=FALSE
                   )
            d
          }
          )

