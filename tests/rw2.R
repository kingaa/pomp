require(pomp)

rw.rprocess <- function (xstart, times, params, ...) { 
  ## this function simulates two independent random walks with intensities s1, s2
  nsims <- ncol(params)
  ntimes <- length(times)
  dt <- diff(times)
  x <- array(0,dim=c(2,nsims,ntimes))
  rownames(x) <- rownames(xstart)
  noise.sds <- params[c('s1','s2'),]
  x[,,1] <- xstart
  for (j in 2:ntimes) {
    ## we are mimicking a continuous-time process, so the increments have SD ~ sqrt(dt)
    ## note that we do not have to assume that 'times' are equally spaced
    x[,,j] <- rnorm(n=2*nsims,mean=x[,,j-1],sd=noise.sds*dt[j-1])
  }
  x
}

rw.dprocess <- function (x, times, params, log = FALSE, ...) { 
  ## given a sequence of consecutive states in 'x', this function computes the p.d.f.
  nsims <- ncol(params)
  ntimes <- length(times)
  dt <- diff(times)
  d <- array(0,dim=c(2,nsims,ntimes-1))
  noise.sds <- params[c('s1','s2'),]
  for (j in 2:ntimes)
    d[,,j-1] <- dnorm(x[,,j]-x[,,j-1],mean=0,sd=noise.sds*dt[j-1],log=TRUE)
  if (log) {
    apply(d,c(2,3),sum)
  } else {
    exp(apply(d,c(2,3),sum))
  }
}


bvnorm.rmeasure <- function (x, times, params, ...) {
  ## noisy observations of the two walks with common noise SD 'tau'
  nsims <- dim(x)[2]
  ntimes <- dim(x)[3]
  y <- array(0,dim=c(2,nsims,ntimes))
  rownames(y) <- c('y1','y2')
  for (k in 1:nsims) {
    for (j in 1:ntimes) {
      y[,k,j] <- rnorm(2,mean=x[,k,j],sd=params['tau',k])
    }
  }
  y
}

bvnorm.dmeasure <- function (y, x, times, params, log = FALSE, ...) {
  ## noisy observations of the two walks with common noise SD 'tau'
  d1 <- dnorm(
              x=y['y1',],
              mean=x['x1',,],
              sd=params['tau',],
              log=TRUE
              )
  d2 <- dnorm(
              x=y['y2',],
              mean=x['x2',,],
              sd=params['tau',],
              log=TRUE
              )
  if (log) {
    d1+d2
  } else {
    exp(d1+d2)
  }
}

rw2 <- pomp(
            rprocess = rw.rprocess,
            dprocess = rw.dprocess,
            rmeasure = bvnorm.rmeasure,
            dmeasure = bvnorm.dmeasure,
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
