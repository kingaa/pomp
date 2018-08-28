library(pomp)
library(magrittr)

pompExample(ou2)

po <- window(ou2,end=10)

set.seed(3434388L)
simulate(po,nsim=5,format="arrays") -> y
y %>% extract2("states") -> x
y %>% extract2("obs") %>% extract(,1,) -> y
t <- time(po)
p <- coef(po)

dmeasure(po,x=x,y=y,times=t,params=p) -> L
dmeasure(po,x=x,y=y,times=t,params=p,log=T) -> ll
stopifnot(
  all.equal(ll[,1:3],log(L[,1:3])),
  identical(dim(ll),c(5L,10L))
)

try(dmeasure("ou2",x=x,y=y,times=t,params=p))
try(dmeasure(x=x,y=y,times=t,params=p))
try(dmeasure(x,y=y,times=t,params=p))
try(dmeasure(po,x=x,y=y,times=t))
try(dmeasure(po,x=x,y=y,params=p))
try(dmeasure(po,x=x,times=t,params=p))
try(dmeasure(po,y=y,times=t,params=p))
try(dmeasure(po,x=as.numeric(x),y=y,times=t,params=p))
try(dmeasure(po,x=x,y=as.numeric(y),times=t,params=p))
try(dmeasure(po,x=x,y=y,times=NULL,params=p))
try(dmeasure(po,x=x[,,1],y=y[,1,drop=FALSE],times=t[1],params=p))
invisible(dmeasure(po,x=x[,,1,drop=FALSE],y=y[,1],times=t[1],params=p))
stopifnot(
  all.equal(dmeasure(po,x=x[,1,,drop=FALSE],y=y,times=t,params=p),
    dmeasure(po,x=x[,1,],y=y,times=t,params=p))
)
try(dmeasure(po,x=x,y=y[1,,drop=FALSE],times=t,params=p))
try(dmeasure(po,x=x[1,,,drop=FALSE],y=y,times=t,params=p))
k <- which(names(p)=="tau")
try(dmeasure(po,x=x,y=y,times=t,params=p[-k]))

pp <- parmat(p,5)
try(dmeasure(po,x=x,y=y,times=t,params=pp[,1:3]))
dmeasure(po,x=x,y=y,times=t,params=pp) -> d
stopifnot(dim(d)==c(5,10),names(dimnames(d))==c("rep","time"))

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

po %>% pomp(
  rmeasure=function(...)c(1,2,3),
  dmeasure=function(...,log)c(3,2)
) -> po1
try(po1 %>% rmeasure(x=x,times=t,params=p))
try(po1 %>% dmeasure(x=x,y=y[,1,],times=t,params=p))

po %>% pomp(
  rmeasure=NULL,
  dmeasure=NULL
) -> po1
po1 %>% rmeasure(x=x,times=t,params=p) %>% is.na() %>% stopifnot()
po1 %>% dmeasure(x=x,y=y[,1,],times=t,params=p) %>% is.na() %>% stopifnot()

pompExample(sir)
po <- window(sir,end=0.5)

set.seed(3434388L)
simulate(po,nsim=5,format="arrays") -> y
y %>% extract2("states") -> x
y %>% extract2("obs") %>% extract(,1,,drop=FALSE) -> y
t <- time(po)
p <- coef(po)

po %>% dmeasure(x=x,y=y,params=p,times=t,log=TRUE) -> d
po %>% rmeasure(x=x,params=p,times=t) -> yy

po %>% pomp(dmeasure=function(...,log)1) %>% dmeasure(x=x,y=y,params=p,times=t) -> d

