options(digits=3)

library(pomp)
library(magrittr)

pompExample(ou2,envir=NULL) -> ou2
ou2[[1]] %>% window(end=10) -> po

set.seed(293982095)

po %>%
  simulate(format="arrays",nsim=3,seed=3434388L) %>%
  extract2("states") -> x
po %>% time() -> t
coef(po) %>% parmat(7) -> p
p["sigma.1",] <- seq(from=1,to=7,by=1)

try(dprocess("ou2",x=x,times=t,params=p))
try(dprocess(x=x,times=t,params=p))
try(po %>% dprocess(x=x,times=t,params=p,log=TRUE))
po %>% dprocess(x=x,times=t,params=p[,1:3],log=TRUE)
try(po %>% dprocess(x=x[,,2],times=t[2],params=p[,1:3],log=FALSE))
try(po %>% dprocess(x=x[,,2:5],times=t[2:5],params=p[,1:2],log=FALSE))
po %>% dprocess(x=x[,,2:5],times=t[2:5],params=p[,1:3],log=TRUE) %>%
  apply(1,sum)
try(po %>% dprocess(x=x[,1:2,2:5],times=t[2:5],params=p[,1:3],log=TRUE) %>%
    apply(1,sum))
po %>% dprocess(x=x[,1,2:5],times=t[2:5],params=p[,1:3],log=TRUE) %>%
  apply(1,sum)

po %>% rinit(params=coef(po)) -> x0
freeze(po %>%
    rprocess(params=coef(po),xstart=parmat(x0,3),times=time(po,t0=TRUE),
      offset=1),
  seed=3434388L) -> x1

stopifnot(max(abs(x-x1))==0)

po %>% rinit(nsim=6) -> x0
try(rprocess("ou2",xstart=x0,times=t,params=p))
try(rprocess(xstart=x0,times=t,params=p))
try(po %>% rprocess(times=t,params=p))
try(po %>% rprocess(xstart=x0,params=p))
try(po %>% rprocess(xstart=x0,times=t))
try(po %>% rprocess(xstart=x0,times=t,params=p))
try(po %>% rprocess(xstart=x0,times=t,params=p[,-1],offset=-3))
try(po %>% rprocess(xstart=x0,times=t,params=p[,-1],offset=500))
try(po %>% rprocess(xstart=x0,times=t,params=p[,-1],offset=NA))
try(po %>% rprocess(xstart=x0,times=t,params=p[,-1],offset=Inf))
po %>% rprocess(xstart=x0,times=t,params=p[,1:3]) -> x
stopifnot(
  dim(x)==c(2,6,10),
  names(dimnames(x))==c("variable","rep","time")
)
po %>% rprocess(xstart=x0[,2],times=t,params=p[,1:3]) -> x
stopifnot(
  dim(x)==c(2,3,10),
  names(dimnames(x))==c("variable","rep","time")
)

try(po %>% rprocess(xstart=x0,times=t[2],params=p))
try(po %>% rprocess(xstart=x0[,2],times=t[2],params=p[,1:3]))
try(po %>% rprocess(xstart=x0[,2:4],times=t[2:5],params=p[,1:2]))
po %>% rprocess(xstart=x0[,2:4],times=t[2:5],params=p[,1:3]) %>%
  apply(1,sum)
