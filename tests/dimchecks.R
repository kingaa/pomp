library(pomp)
library(reshape2)
library(plyr)
library(magrittr)

set.seed(1420306530L)

pompExample(ricker)
po <- ricker

pars <- coef(po)
xstart <- init.state(po,params=pars)
rprocess(po,xstart,times=0:5,params=pars)[,1,]
rprocess(po,xstart=parmat(xstart,5),times=0:5,params=pars)[,3,]
rprocess(po,xstart=xstart,times=0:5,params=parmat(pars,3))[,3,]
try(rprocess(po,xstart=parmat(xstart,2),times=0:5,params=parmat(pars,3))[,,3])
rprocess(po,xstart=parmat(xstart,2),times=0:5,params=parmat(pars,6))[,,3]
x <- rprocess(po,xstart=parmat(xstart,2),times=0:5,params=parmat(pars,8))
rmeasure(po,x=x,params=pars,times=0:5)[,,3]
try(rmeasure(po,x=x,params=parmat(pars,3),times=0:5)[,,3])
rmeasure(po,x=x,params=parmat(pars,4),times=0:5)
x <- rprocess(po,xstart=xstart,times=0:5,params=pars)
rmeasure(po,x=x,params=parmat(pars,2),times=0:5)
y <- rmeasure(po,x=x,params=parmat(pars,4),times=0:5)
try(rmeasure(po,x=x,params=parmat(pars,4),times=numeric(0)))
try(rmeasure(po,x=x,params=parmat(pars,4),times=0:2))
dmeasure(po,x=x,y=y[,2,,drop=F],params=pars,times=0:5)
x <- rprocess(po,xstart=parmat(xstart,3),times=0:5,params=pars)
y <- rmeasure(po,x=x,params=pars,times=0:5)
try(dmeasure(po,x=x,y=y,params=parmat(pars,3),times=0:5))
try(dmeasure(po,x=x[,,1:3,drop=F],y=y[,2,,drop=F],params=parmat(pars,3),times=0:5))
try(dmeasure(po,x=x,y=y[,1,,drop=F],params=parmat(pars,2),times=0:5))
try(dmeasure(po,x=x[,1:2,],y=y[,1,,drop=F],params=parmat(pars,3),times=0:5))
f1 <- dmeasure(po,x=x,y=y[,1,,drop=F],params=parmat(pars,3),times=0:5)
f2 <- dmeasure(po,x=x,y=y[,1,,drop=F],params=pars,times=0:5)
stopifnot(identical(f1,f2))
try(dmeasure(po,x=x,y=y,params=parmat(pars,3),times=numeric(0)))

g1 <- skeleton(po,x=x,t=0:5,params=pars)
g2 <- skeleton(po,x=x,t=0:5,params=parmat(pars,3))
stopifnot(identical(g1,g2))
g3 <- skeleton(po,x=x,t=0:5,params=parmat(pars,6))
stopifnot(identical(g1,g3[,1:3,]))
stopifnot(identical(g1,g3[,4:6,]))

pompExample(gompertz)
p <- parmat(coef(gompertz),5)
f1 <- partrans(gompertz,p,"inv")
f2 <- parmat(coef(gompertz,transform=TRUE),5)
stopifnot(identical(f1,f2))

pars1 <- parmat(coef(ricker),3)
pars1["N.0",] <- c(3,5,7)
try(xstart <- init.state(po,params=pars1,nsim=8))
xstart <- init.state(po,params=pars1,nsim=6)
rprocess(po,xstart,times=0:5,params=pars1)[1,4:6,6]

pars <- coef(ricker)
pars["r"] <- 10
xstart <- init.state(ricker,params=pars,nsim=8)
rprocess(po,xstart,times=0:5,params=pars)[,1,]
simulate(ricker,params=pars1,nsim=2,times=0) %>%
    ldply(as.data.frame)
simulate(ricker,params=pars1,nsim=1,times=0:1,as=T) %>%
    mutate(N=signif(N,3),e=signif(e,3))
simulate(ricker,params=pars1,nsim=1,as=T,include.data=T) %>%
    melt(id=c("time","sim")) %>%
    na.omit() %>%
    ddply(~variable+sim,summarize,mean=signif(mean(value),2))
simulate(ricker,params=pars1,nsim=2,times=0:1,states=T) %>%
    melt() %>% dcast(rep+time~variable) %>%
    mutate(N=signif(N,3),e=signif(e,3))
simulate(ricker,params=pars1,nsim=2,times=0:1,obs=T) %>%
    melt() %>% dcast(rep+time~variable)
simulate(ricker,params=pars1,nsim=2,times=0:1,obs=T,states=T) %>%
    melt() %>% dcast(rep+time~variable) %>%
    mutate(N=signif(N,3),e=signif(e,3))

pomp(ricker,
     initializer=Csnippet("N = rnorm(7,1); e = 0;"),statenames=c("N","e")
     ) -> ricker
xstart <- init.state(ricker,nsim=6)
pars1 <- parmat(coef(ricker),3)
pars1["r",] <- c(2,5,10)
pars1["sigma",] <- 0.0
x <- rprocess(ricker,params=pars1,times=c(0,1),xstart=xstart)
stopifnot(all.equal(pars1["r",]*xstart[1,]*exp(-xstart[1,]),x["N",,2]))

pompExample(ou2)

stopifnot(identical(dim(x0 <- init.state(ou2,params=coef(ou2))),c(2L,1L)))
stopifnot(identical(dim(init.state(ou2,params=parmat(coef(ou2),5))),c(2L,5L)))
try(rprocess(ou2,xstart=x0,times=0,params=coef(ou2),offset=0))
try(rprocess(ou2,xstart=x0,times=seq(0,10),params=coef(ou2),offset=-5))
stopifnot(identical(
    dim(rprocess(ou2,xstart=x0,times=seq(0,10),params=coef(ou2),offset=10)),
    c(2L,1L,1L)))
try(rprocess(ou2,xstart=x0,times=seq(0,10),params=coef(ou2),offset=11))
stopifnot(identical(
    dim(rprocess(ou2,xstart=x0,times=seq(0,10),params=coef(ou2),offset=0,ignored=3)),
    c(2L,1L,11L)))
stopifnot(identical(
    dim(rprocess(ou2,xstart=parmat(x0,3),times=seq(0,10),params=coef(ou2),offset=1)),
    c(2L,3L,10L)))
stopifnot(identical(
    dim(rprocess(ou2,xstart=parmat(x0,3),times=seq(0,10),params=parmat(coef(ou2),3),offset=1)),
    c(2L,3L,10L)))
stopifnot(identical(
    dim(x <- rprocess(ou2,xstart=parmat(x0,3),times=seq(0,10),params=parmat(coef(ou2),6),offset=1)),
    c(2L,6L,10L)))
try(rprocess(ou2,xstart=parmat(x0,3),times=seq(0,10),params=parmat(coef(ou2),5),offset=1))
try(rprocess(ou2,xstart=parmat(x0,11),times=seq(0,10),params=parmat(coef(ou2),5),offset=1))

pomp(ou2,rprocess=function(xstart,times,params,...){xstart}) -> po
try(rprocess(po,xstart=parmat(x0,3),times=seq(0,10),params=parmat(coef(ou2),1),offset=1))
pomp(ou2,rprocess=function(xstart,times,params,...){x <- xstart; dim(x) <- c(nrow(xstart),1,ncol(xstart)); x}) -> po
try(rprocess(po,xstart=parmat(x0,3),times=seq(0,10),params=parmat(coef(ou2),1),offset=1))
pomp(ou2,rprocess=function(xstart,times,params,...){array(runif(66),dim=c(2,3,11))}) -> po
try(rprocess(po,xstart=parmat(x0,3),times=seq(0,10),params=parmat(coef(ou2),1),offset=1))

try(dprocess(ou2,x=x[,3,1],times=1,params=coef(ou2),log=TRUE))
stopifnot(identical(
    dim(dprocess(ou2,x=x[,3,1:3],times=1:3,params=coef(ou2),log=TRUE)),
    c(1L,2L)))
stopifnot(identical(
    dim(dprocess(ou2,x=x[,3:5,],times=1:10,params=coef(ou2),log=TRUE)),
    c(3L,9L)))
try(dprocess(ou2,x=x[,3:5,],times=1:8,params=coef(ou2),log=TRUE))
stopifnot(identical(
    dim(dprocess(ou2,x=x[,3:5,],times=1:10,params=parmat(coef(ou2),1))),
    c(3L,9L)))
stopifnot(identical(
    dim(dprocess(ou2,x=x[,3:5,],times=1:10,params=parmat(coef(ou2),3))),
    c(3L,9L)))
try(dprocess(ou2,x=x[,3:5,],times=1:10,params=parmat(coef(ou2),2)))
stopifnot(identical(
    dim(dprocess(ou2,x=x[,3:6,],times=1:10,params=parmat(coef(ou2),2))),
    c(4L,9L)))
try(dprocess(ou2,x=x[,3:6,],times=1:10,params=parmat(coef(ou2),3)))
stopifnot(identical(
    dim(dprocess(ou2,x=x[,3:6,],times=1:10,params=parmat(coef(ou2),4))),
    c(4L,9L)))

pomp(ou2,dprocess=function(x,times,params,log,...){0}) -> po
try(dprocess(po,x=x[,3:6,],times=1:10,params=parmat(coef(ou2),4)))
pomp(ou2,dprocess=function(x,times,params,log,...){array(0,dim=c(1,1))}) -> po
try(dprocess(po,x=x[,3:6,],times=1:10,params=parmat(coef(ou2),4)))
pomp(ou2,dprocess=function(x,times,params,log,...){array(0,dim=c(dim(x)[2],length(times)-1))}) -> po
stopifnot(identical(
    dim(dprocess(po,x=x[,3:6,],times=1:10,params=parmat(coef(ou2),4))),
    c(4L,9L)))
stopifnot(identical(
    dim(dprocess(po,x=x,times=1:10,params=parmat(coef(ou2),2))),
    c(6L,9L)))

try(simulate(ou2,nsim=0))
try(class(simulate(ou2,nsim=c(1,20))))
try(simulate(ou2,nsim=numeric(0)))

po <- simulate(ou2,times=1:5)
xs <- xstart <- init.state(po)
tt <- t <- time(po)
xx <- x <- states(po)
pp <- p <- coef(po)
yy <- y <- obs(po)
storage.mode(xs) <- "integer"
storage.mode(xx) <- "integer"
storage.mode(pp) <- "integer"
storage.mode(tt) <- "integer"
storage.mode(yy) <- "integer"
dim(rprocess(po,xstart=xs,times=t,params=p))
dim(rprocess(po,xstart=xstart,times=t,params=pp))
dim(rprocess(po,xstart=xs,times=t,params=p))
dim(rprocess(po,xstart=xstart,times=tt,params=p))
dim(rmeasure(po,x=xx,params=p,times=t))
dim(rmeasure(po,x=x,params=pp,times=t))
dim(rmeasure(po,x=x,params=p,times=tt))
dim(dprocess(po,x=xx,params=p,times=t))
dim(dprocess(po,x=x,params=pp,times=t))
dim(dprocess(po,x=x,params=p,times=tt))
dim(dmeasure(po,y=yy,x=x,params=p,times=t))
dim(dmeasure(po,y=y,x=xx,params=p,times=t))
dim(dmeasure(po,y=y,x=x,params=pp,times=t))
dim(dmeasure(po,y=y,x=x,params=p,times=tt))
dim(init.state(po,params=p,t0=t))
dim(init.state(po,params=pp,t0=t))
dim(init.state(po,params=p,t0=tt))
dim(skeleton(po,x=xx,t=t,params=p))
dim(skeleton(po,x=x,t=tt,params=p))
dim(skeleton(po,x=x,t=t,params=pp))
dim(trajectory(po,times=tt,params=p))
dim(trajectory(po,times=t,params=pp))
length(dprior(po,params=p))
length(dprior(po,params=pp))

