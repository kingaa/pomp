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

bad.initializer <- function (params, t0, ...) 
{
    ivpnames <- c("x1.0","x2.0")
    x <- params[ivpnames]
    x
}

crap.initializer <- function (params, t0, ...) 
{
    x <- rnorm(n=ceiling(runif(n=1,min=0,max=10)))
    names(x) <-head(letters,length(x))
    x
}

p <- rbind(s1=c(2,2,3),s2=c(0.1,1,2),tau=c(1,5,0),x1.0=c(0,0,5),x2.0=c(0,0,0))

rw2 <- pomp(
    rprocess = rw.rprocess,
    dprocess = rw.dprocess,
    measurement.model=list(
        y1 ~ norm(mean=x1,sd=tau),
        y2 ~ norm(mean=x2,sd=tau)
    ),
    
    initializer=bad.initializer,
    times=1:100,
    data=rbind(
        y1=rep(0,100),
        y2=rep(0,100)
    ),
    t0=0,
    useless=23
)

invisible(show(rw2))

try(
    simulate(rw2,params=p)
)

rw2 <- pomp(rw2,initializer=crap.initializer)

try(
    simulate(rw2,params=p)
)

rw2 <- pomp(
    rprocess = rw.rprocess,
    dprocess = rw.dprocess,
    measurement.model=list(
        y1 ~ norm(mean=x1,sd=tau),
        y2 ~ norm(mean=x2,sd=tau)
    ),
    rmeasure=Csnippet("sid"),
    dmeasure=Csnippet("nancy"),
    times=1:100,
    data=rbind(
        y1=rep(0,100),
        y2=rep(0,100)
    ),
    t0=0
)

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

x <- obs(rw2)
t <- time(rw2)

x0 <- init.state(rw2,params=p)
x <- rprocess(rw2,xstart=x0,times=0:100,params=p)
y <- rmeasure(rw2,x=x,times=0:100,params=p)

a1 <- dmeasure(rw2,y=y[,1,1:4],x=x[,,1:4,drop=F],times=time(rw2)[1:4],p)
b1 <- dmeasure(rw2,y=y[,2,1:4],x=x[,,1:4,drop=F],times=time(rw2)[1:4],p)
c1 <- log(dmeasure(rw2,y=y[,3,1:4],x=x[,,1:4,drop=F],times=time(rw2)[1:4],p))
d1 <- dmeasure(rw2,y=y[,3,1:4],x=x[,,1:4,drop=F],times=time(rw2)[1:4],p,log=T)
e1 <- dprocess(rw2,x[,,6:11],times=5:10,params=p,log=T)
f1 <- log(dprocess(rw2,x[,,6:11],times=5:10,params=p))
stopifnot(max(abs(c1-d1),na.rm=T)<.Machine$double.eps*100)
stopifnot(max(abs(e1-f1),na.rm=T)<.Machine$double.eps*100)

po <- pomp(
    rw2,
    dmeasure = bvnorm.dmeasure,
    rmeasure = bvnorm.rmeasure
)

a2 <- dmeasure(po,y=y[,1,1:4],x=x[,,1:4,drop=F],times=time(rw2)[1:4],p)
b2 <- dmeasure(po,y=y[,2,1:4],x=x[,,1:4,drop=F],times=time(rw2)[1:4],p)
c2 <- log(dmeasure(po,y=y[,3,1:4],x=x[,,1:4,drop=F],times=time(rw2)[1:4],p))
d2 <- dmeasure(po,y=y[,3,1:4],x=x[,,1:4,drop=F],times=time(rw2)[1:4],p,log=T)
e2 <- dprocess(po,x[,,6:11],times=5:10,params=p,log=T)
f2 <- log(dprocess(rw2,x[,,6:11],times=5:10,params=p))
stopifnot(max(abs(c2-d2),na.rm=T)<.Machine$double.eps*100)
stopifnot(max(abs(e2-f2),na.rm=T)<.Machine$double.eps*100)

stopifnot(max(abs(a1-a2),na.rm=T)<.Machine$double.eps*100)
stopifnot(max(abs(b1-b2),na.rm=T)<.Machine$double.eps*100)
stopifnot(max(abs(d1-d2),na.rm=T)<.Machine$double.eps*100)
stopifnot(max(abs(e1-e2),na.rm=T)<.Machine$double.eps*100)

pompExample(rw2)

a3 <- dmeasure(po,y=y[,1,1:4],x=x[,,1:4,drop=F],times=time(rw2)[1:4],p)
b3 <- dmeasure(po,y=y[,2,1:4],x=x[,,1:4,drop=F],times=time(rw2)[1:4],p)
c3 <- log(dmeasure(po,y=y[,3,1:4],x=x[,,1:4,drop=F],times=time(rw2)[1:4],p))
d3 <- dmeasure(po,y=y[,3,1:4],x=x[,,1:4,drop=F],times=time(rw2)[1:4],p,log=T)
e3 <- dprocess(po,x[,,6:11],times=5:10,params=p,log=T)
f3 <- log(dprocess(rw2,x[,,6:11],times=5:10,params=p))
stopifnot(max(abs(c3-d3),na.rm=T)<.Machine$double.eps*100)
stopifnot(max(abs(e3-f3),na.rm=T)<.Machine$double.eps*100)

stopifnot(max(abs(a2-a3),na.rm=T)<.Machine$double.eps*100)
stopifnot(max(abs(b2-b3),na.rm=T)<.Machine$double.eps*100)
stopifnot(max(abs(d2-d3),na.rm=t)<.Machine$double.eps*100)
stopifnot(max(abs(e2-e3),na.rm=t)<.Machine$double.eps*100)

new <- window(rw2,start=20,end=30)
new <- simulate(new)

invisible(timezero(new))
timezero(new) <- 19
print(simulate(new))

time(rw2) <- seq(1,1000,by=20)
x <- simulate(rw2)
invisible(states(x)[,1:5])
try(
    time(rw2) <- seq(-20,1000,by=20)
)
try(
    time(rw2) <- c(0,5,10,15,12,20)
)
time(rw2,t0=TRUE) <- seq(-20,1000,by=20)
x <- simulate(rw2)
time(rw2) <- c(0,20,25.8,50,60)
time(rw2,t0=TRUE) <- c(0,20,25.8,50,60)
time(rw2,t0=TRUE) <- c(0,0,20,25.8,50,60)
time(rw2) <- c(0,20,25.8,50,60)
