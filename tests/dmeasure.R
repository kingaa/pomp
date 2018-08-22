library(pomp)
library(magrittr)
library(plyr)
library(reshape2)

pompExample(ou2)

po <- window(ou2,end=10)

set.seed(3434388L)
simulate(po,nsim=5,as.data.frame=TRUE,include.data=FALSE)%>%
  melt(id.vars=c("time","sim")) %>%
  acast(variable~sim~time) -> y
x <- y[c("x1","x2"),,,drop=FALSE]
y <- y[c("y1","y2"),1,]
t <- time(po)
theta <- coef(po)

dmeasure(po,x=x,y=y,times=t,params=theta) -> L
dmeasure(po,x=x,y=y,times=t,params=theta,log=T) -> ll
stopifnot(all.equal(ll[,1:3],log(L[,1:3])))
stopifnot(identical(dim(ll),c(5L,10L)))

try(dmeasure(x=x,y=y,times=t,params=theta))
try(dmeasure(x,y=y,times=t,params=theta))
try(dmeasure(po,x=x,y=y,times=t))
try(dmeasure(po,x=x,y=y,params=theta))
try(dmeasure(po,x=x,times=t,params=theta))
try(dmeasure(po,y=y,times=t,params=theta))
try(dmeasure(po,x=as.numeric(x),y=y,times=t,params=theta))
try(dmeasure(po,x=x,y=as.numeric(y),times=t,params=theta))
try(dmeasure(po,x=x,y=y,times=NULL,params=theta))
try(dmeasure(po,x=x[,,1],y=y[,1,drop=FALSE],times=t[1],params=theta))
invisible(dmeasure(po,x=x[,,1,drop=FALSE],y=y[,1],times=t[1],params=theta))
stopifnot(all.equal(dmeasure(po,x=x[,1,,drop=FALSE],y=y,times=t,params=theta),
                    dmeasure(po,x=x[,1,],y=y,times=t,params=theta)))
try(dmeasure(po,x=x,y=y[1,,drop=FALSE],times=t,params=theta))
try(dmeasure(po,x=x[1,,,drop=FALSE],y=y,times=t,params=theta))
k <- which(names(theta)=="tau")
try(dmeasure(po,x=x,y=y,times=t,params=theta[-k]))

po %>% pomp(dmeasure=NULL) -> po1
stopifnot(sum(is.na(dmeasure(po1,x=x,y=y,times=t,params=theta)))==50)

xx <- x
xx[] <- as.integer(round(x))
stopifnot(all.equal(dmeasure(po,x=xx,y=y,times=t,params=theta,log=TRUE),
                    dmeasure(po,x=round(x),y=y,times=t,params=theta,log=TRUE)))

stopifnot({
  dmeasure(pomp(po,dmeasure=function(y,x,t,params,log,...)3),
           x=x,y=y,times=t,params=theta,log=TRUE) -> z
  dimnames(z)
  all(z[]==3)
})

try(dmeasure(pomp(po,dmeasure=function(y,x,t,params,log,...)c(3,2,1)),
             x=x,y=y,times=t,params=theta,log=TRUE))

pp <- parmat(theta,10)
dmeasure(po,params=pp,x=x,y=y,t=t,log=TRUE) -> ll
stopifnot(all.equal(ll[1:5,],ll[6:10,]))
try(dmeasure(po,params=pp[,1:7],x=x,y=y,t=t,log=TRUE))

pompExample(dacca)
set.seed(3434388L)
po <- window(dacca,end=1892)
po %>% simulate(nsim=5) -> dat
dat %>% lapply(states) %>% melt() %>%
  acast(variable~L1~time) -> x
dat %>% lapply(obs) %>% melt() %>%
  subset(L1==1) %>%
  acast(variable~L1~time) -> y
t <- time(po)
theta <- coef(po)
dmeasure(po,x=x,y=y,times=t,params=theta) -> L
dmeasure(po,x=x,y=y,times=t,params=theta,log=T) -> ll
stopifnot(all.equal(ll[,1:3],log(L[,1:3])))
