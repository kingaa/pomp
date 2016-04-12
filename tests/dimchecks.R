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

try(
    rprocess(po,xstart=parmat(xstart,2),times=0:5,params=parmat(pars,3))[,,3]
    )

rprocess(po,xstart=parmat(xstart,2),times=0:5,params=parmat(pars,6))[,,3]

x <- rprocess(po,xstart=parmat(xstart,2),times=0:5,params=parmat(pars,8))

rmeasure(po,x=x,params=pars,times=0:5)[,,3]

try(
    rmeasure(po,x=x,params=parmat(pars,3),times=0:5)[,,3]
    )

rmeasure(po,x=x,params=parmat(pars,4),times=0:5)

x <- rprocess(po,xstart=xstart,times=0:5,params=pars)
rmeasure(po,x=x,params=parmat(pars,2),times=0:5)

y <- rmeasure(po,x=x,params=parmat(pars,4),times=0:5)
dmeasure(po,x=x,y=y[,2,,drop=F],params=pars,times=0:5)

x <- rprocess(po,xstart=parmat(xstart,3),times=0:5,params=pars)
y <- rmeasure(po,x=x,params=pars,times=0:5)
try(dmeasure(po,x=x,y=y,params=parmat(pars,3),times=0:5))
f1 <- dmeasure(po,x=x,y=y[,1,,drop=F],params=parmat(pars,3),times=0:5)
f2 <- dmeasure(po,x=x,y=y[,1,,drop=F],params=pars,times=0:5)
stopifnot(identical(f1,f2))

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
