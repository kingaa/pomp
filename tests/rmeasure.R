options(digits=3)

library(pomp)

pompExample(ou2)

po <- window(ou2,end=10)

set.seed(3434388L)
x <- simulate(po,obs=F,states=T,nsim=100)
t <- time(po)
theta <- coef(po)

rmeasure(po,x=x,times=t,params=theta) -> y
stopifnot(identical(dim(y),c(2L,100L,10L)))
fivenum(as.numeric(y-x)/coef(po,"tau"))

try(rmeasure(x=x,times=t,params=theta))
try(rmeasure(po,x=x,times=t))
try(rmeasure(po,x=x,params=theta))
try(rmeasure(po,times=t,params=theta))
try(rmeasure(po,x=as.numeric(x),y=y,times=t,params=theta))
try(rmeasure(po,x=x,times=NULL,params=theta))
try(rmeasure(po,x=x[,,1],times=t[1],params=theta))
invisible(rmeasure(po,x=x[,,1,drop=FALSE],times=t[1],params=theta))
freeze(rmeasure(po,x=x[,1,],times=t,params=theta),seed=349596L) -> y1
freeze(rmeasure(po,x=x[,1,,drop=FALSE],times=t,params=theta),seed=349596L) -> y2
stopifnot(all.equal(as.numeric(y1),as.numeric(y2)))
try(rmeasure(po,x=x[1,,,drop=FALSE],times=t,params=theta))
k <- which(names(theta)=="tau")
try(rmeasure(po,x=x,y=y,times=t,params=theta[-k]))

stopifnot({
  rmeasure(pomp(po,rmeasure=function(y,x,t,params,log,...)c(3,3)),
    x=x,y=y,times=t,params=theta,log=TRUE) -> z
  dimnames(z)
  all(z[]==3)
})

try(rmeasure(pomp(po,rmeasure=function(x,t,params,...)c(3,2,1)),
  x=x,y=y,times=t,params=theta,log=TRUE))

pp <- parmat(theta,10)
rmeasure(po,params=pp,x=x,y=y,t=t,log=TRUE) -> y
try(rmeasure(po,params=pp[,1:7],x=x,t=t,log=TRUE))

pompExample(dacca)
set.seed(3434388L)
po <- window(dacca,end=1892)
dat <- simulate(po,obs=T,states=T,nsim=2)
