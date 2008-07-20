library(pomp)

rw.rprocess <- function (params, xstart, times, ...) { 
  ## this function simulates two independent random walks with intensities s1, s2
  nvars <- nrow(xstart)
  nreps <- ncol(params)
  ntimes <- length(times)
  dt <- diff(times)
  x <- array(0,dim=c(nvars,nreps,ntimes))
  rownames(x) <- rownames(xstart)
  noise.sds <- params[c('s1','s2'),]
  x[,,1] <- xstart
  for (j in 2:ntimes) {
    ## we are mimicking a continuous-time process, so the increments have SD ~ sqrt(dt)
    ## note that we do not have to assume that 'times' are equally spaced
    x[c("x1","x2"),,j] <- rnorm(
                                n=2*nreps,
                                mean=x[c("x1","x2"),,j-1],
                                sd=noise.sds*dt[j-1]
                                )
  }
  x
}

rw.dprocess <- function (x, times, params, log = FALSE, ...) { 
  ## given a sequence of consecutive states in 'x', this function computes the p.d.f.
  nreps <- ncol(params)
  ntimes <- length(times)
  dt <- diff(times)
  d <- array(0,dim=c(2,nreps,ntimes-1))
  noise.sds <- params[c('s1','s2'),]
  for (j in 2:ntimes)
    d[,,j-1] <- dnorm(x[,,j]-x[,,j-1],mean=0,sd=noise.sds*dt[j-1],log=TRUE)
  d <- apply(d,c(2,3),sum)
  if (log) d else exp(d)
}

bvnorm.rmeasure <- function (t, x, params, ...) {
  ## noisy observations of the two walks with common noise SD 'tau'
  c(
    y1=rnorm(n=1,mean=x['x1'],sd=params['tau']),
    y2=rnorm(n=1,mean=x['x2'],sd=params['tau'])
    )
}

bvnorm.dmeasure <- function (y, x, t, params, log = FALSE, ...) {
  f <- sum(
           dnorm(
                 x=y[c("y1","y2")],
                 mean=x[c("x1","x2")],
                 sd=params["tau"],
                 log=TRUE
                 ),
           na.rm=TRUE
           )
  if (log) f else exp(f)
}

rw2 <- pomp(
            rprocess = rw.rprocess,
            dprocess = rw.dprocess,
            measurement.model=list(
              y1 ~ norm(mean=x1,sd=tau),
              y2 ~ norm(mean=x2,sd=tau)
            ),
            times=1:100,
            data=rbind(
              y1=rep(0,100),
              y2=rep(0,100)
              ),
            t0=0,
            useless=23
            )

po <- pomp(
           rprocess = euler.simulate,
           dprocess = euler.density,
           delta.t = 1,
           step.fun = function(x, t, params, dt, ...) {
             c(
               y1=rnorm(n=1,mean=x['x1'],sd=params['s1']),
               y2=rnorm(n=1,mean=x['x2'],sd=params['s2'])
               )
           },
           dens.fun = function (x1, t1, x2, t2, params, ...) {
             sum(
                 dnorm(
                       x=x2[c('x1','x2')],
                       mean=x1[c('x1','x2')],
                       sd=params[c('s1','s2')]
                       ),
                 na.rm=TRUE
                 )
           },
           measurement.model=list(
             y1 ~ norm(mean=x1,sd=tau),
             y2 ~ norm(mean=x2,sd=tau)
             ),
           times=1:100,
           data=rbind(
             y1=rep(0,100),
             y2=rep(0,100)
             ),
           t0=0
           )
