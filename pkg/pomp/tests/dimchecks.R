library(pomp)

set.seed(1420306530L)

data(ricker)
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

data(gompertz)
p <- parmat(coef(gompertz),5)
f1 <- partrans(gompertz,p,"inv")
f2 <- parmat(coef(gompertz,transform=TRUE),5)
stopifnot(identical(f1,f2))

