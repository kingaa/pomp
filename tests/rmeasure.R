library(pomp)
library(magrittr)

pompExample(ou2)

po <- window(ou2,end=10)

set.seed(3434388L)
simulate(po,nsim=5,format="arrays") -> y
y %>% extract2("states") -> x
t <- time(po)
p <- coef(po)

rmeasure(po,x=x,times=t,params=p) -> y
stopifnot(
  dim(y)==c(2,5,10),
  names(dimnames(y))==c("variable","rep","time")
)

try(rmeasure("ou2",x=x,times=t,params=p))
try(rmeasure(x=x,times=t,params=p))
try(rmeasure(x,times=t,params=p))
try(rmeasure(po,x=x,times=t))
try(rmeasure(po,x=x,params=p))
try(rmeasure(po,x=x,times=t,params=p))
try(rmeasure(po,x=as.numeric(x),times=t,params=p))
try(rmeasure(po,x=x,times=NULL,params=p))
try(rmeasure(po,x=x[,,1],times=t[1],params=p))
invisible(rmeasure(po,x=x[,,1,drop=FALSE],times=t[1],params=p))
try(rmeasure(po,x=x[1,,,drop=FALSE],times=t,params=p))
k <- which(names(p)=="tau")
try(rmeasure(po,x=x,y=y,times=t,params=p[-k]))

pp <- parmat(p,5)
try(rmeasure(po,x=x,times=t,params=pp[,1:3]))
rmeasure(po,x=x,times=t,params=pp) -> y
stopifnot(dim(y)==c(2,5,10),names(dimnames(y))==c("variable","rep","time"))
