library(pomp)

set.seed(45768683)

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

p <- rbind(s1=c(2,2,3),s2=c(0.1,1,2),tau=c(1,5,0),x1.0=c(0,0,5),x2.0=c(0,0,0))
examples <- simulate(rw2,params=p)
rw2 <- examples[[1]]

y <- simulate(rw2,params=p,obs=T,states=T)
y <- simulate(rw2,params=p,obs=T)
x <- simulate(rw2,params=p,states=T)
x <- simulate(rw2,nsim=10,params=p,states=T)
x <- simulate(rw2,nsim=10,params=p[,1],states=T)
x <- simulate(rw2,nsim=10,params=p[,1],obs=T,states=T)
x <- simulate(rw2,nsim=10,params=p[,1],obs=T,states=T)
x <- simulate(rw2,nsim=10,params=p,obs=T,states=T)
x <- simulate(rw2,nsim=10,params=p[,1])

x <- data.array(rw2)
t <- time(rw2)

x0 <- init.state(rw2,params=p)
x <- rprocess(rw2,xstart=x0,times=0:100,params=p)
y <- rmeasure(rw2,x=x,times=0:100,params=p)

log(dprocess(rw2,x[,,6:11],times=5:10,params=p))
dprocess(rw2,x[,,6:11],times=5:10,params=p,log=T)

dmeasure(rw2,y=y[,1,1:4],x=x[,,1:4,drop=F],times=time(rw2)[1:4],p)
dmeasure(rw2,y=y[,2,1:4],x=x[,,1:4,drop=F],times=time(rw2)[1:4],p)
log(dmeasure(rw2,y=y[,3,1:4],x=x[,,1:4,drop=F],times=time(rw2)[1:4],p))
dmeasure(rw2,y=y[,3,1:4],x=x[,,1:4,drop=F],times=time(rw2)[1:4],p,log=T)

po <- pomp(
           rprocess = rw.rprocess,
           dprocess = rw.dprocess,
           dmeasure = bvnorm.dmeasure,
           rmeasure = bvnorm.rmeasure,
           times=1:100,
           data=rbind(
             y1=rep(0,100),
             y2=rep(0,100)
             ),
           t0=0
           )

dmeasure(po,y=y[,1,1:4],x=x[,,1:4,drop=F],times=time(rw2)[1:4],p)
dmeasure(po,y=y[,2,1:4],x=x[,,1:4,drop=F],times=time(rw2)[1:4],p)
log(dmeasure(po,y=y[,3,1:4],x=x[,,1:4,drop=F],times=time(rw2)[1:4],p))
dmeasure(po,y=y[,3,1:4],x=x[,,1:4,drop=F],times=time(rw2)[1:4],p,log=T)

data(rw2)

dmeasure(po,y=y[,1,1:4],x=x[,,1:4,drop=F],times=time(rw2)[1:4],p)
dmeasure(po,y=y[,2,1:4],x=x[,,1:4,drop=F],times=time(rw2)[1:4],p)
dmeasure(po,y=y[,3,1:4],x=x[,,1:4,drop=F],times=time(rw2)[1:4],p,log=T)
dprocess(rw2,x[,,6:11],times=5:10,params=p,log=T)

