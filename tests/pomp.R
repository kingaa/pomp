options(digits=3)

set.seed(758723694)

library(pomp)
library(magrittr)

try(pomp())
try(pomp("bob"))
try(pomp(times=3))
stopifnot(is.null(pomp(NULL)))
try(data.frame(a=1:10,a=1:10,check.names=FALSE) %>% pomp(t0=4))
try(data.frame(a=1:10,b=1:10) %>% pomp(t0=4))
try(data.frame(a=1:10,b=1:10) %>% pomp(times="b"))
try(data.frame(a=10:1,b=1:10) %>% pomp(times="a",t0=0))
try(data.frame(a=1:10,b=1:10) %>% pomp(times=1:10))
try(data.frame(a=1:10,b=1:10) %>% pomp(times=1))
try(data.frame(a=1:10,b=1:10) %>% pomp(times=3))
try(data.frame(a=1:10,b=1:10) %>% pomp(times="c"))
try(data.frame(a=1:10,b=1:10) %>% pomp(times=NA))
try(data.frame(a=1:10,b=1:10) %>% pomp(times=NULL))
try(data.frame(a=1:10,b=1:10) %>% pomp(times=1,t0=11))
try(data.frame(a=1:10,b=1:10) %>% pomp(times=1,t0=NULL))
try(data.frame(a=1:10,b=1:10) %>% pomp(times=1,t0=NA))
stopifnot(data.frame(a=1:10,b=1:10) %>%
    pomp(covar=covariate_table(c=0:10,d=0:10,times=1),
      covarnames="d",times=1,t0=0,bob=3) %>% class() %>%
    equals("pomp"))
try(data.frame(a=1:10,b=1:10) %>%
    pomp(covar=covariate_table(c=1:10,d=1:10,d=1:10,times=1),
      times=1,t0=0))
stopifnot(data.frame(a=1:10,b=1:10) %>%
    pomp(times="a",t0=0) %>% class() %>%
    equals("pomp"))
try(simulate())
try(NULL %>% pomp(t0=4))
try(NULL %>% pomp(times="a",t0=0))
try(NULL %>% pomp(times=1:10,t0=3))
try(NULL %>% pomp(times=1:10,t0=1,rinit=3))
stopifnot(NULL %>% pomp(times=1:10,t0=1) %>% class() %>% equals("pomp"))

pompExample(gompertz,envir=NULL) %>% extract2(1) -> po
try(po %>% pomp(rprocess=NULL) %>% simulate())
po %>% pomp(rprocess=NULL) %>% slot("rprocess")
po %>% pomp(skeleton=NULL) %>% slot("skeleton")
po %>% pomp(partrans=NULL) %>% slot("partrans")
stopifnot({
  po %>% pomp(partrans=NULL) %>% coef(transform=TRUE) -> theta1
  coef(po) -> theta2
  theta1==theta2
},
  po %>% pomp(times=1:5) %>% class() %>% equals("pomp"))

stopifnot(po %>%
  pomp(rprocess=onestep.sim(function(x,t,params,delta.t,...)x),
    skeleton=map(function(x,t,params,...)x),
    rmeasure=function(x,t,params,...)3,
    dmeasure=function(x,y,t,params,log,...)1,
    covar=covariate_table(a=1:20,b=1:20,times="a")) %>% class() %>%
  equals("pomp"))

try(po %>% pomp(times=3:1))

try(po %>% pomp(rinit=Csnippet("X=3;")))
stopifnot(
  po %>% pomp(rinit=Csnippet("X=3;"),statenames=c("X","Z")) %>%
    class() %>% equals("pomp"))
try(po %>% pomp(rprocess="bob"))
try(po %>% pomp(skeleton="bob"))
try(po %>% pomp(partrans="bob"))
try(po %>% pomp(params=c(1,2,3)))
try(po %>% pomp(params=c(a=1,b=2,3)))

coef(po)
coef(po,transform=TRUE)
coef(po,c("r","tau"))
try(coef(po,c("bob","tau")))
try(coef(po) <- c(1,2,3))
try(coef(po,transform=TRUE) <- c(1,2,3))
coef(po) <- list(as.list(coef(po)))
coef(po,"r") <- 0.2
coef(po,"r") <- list(r=0.2)
coef(po,c("r","theta")) <- list(r=0.2)
coef(po,"sigma",transform=TRUE) <- 0
coef(po)
coef(po) <- NULL
stopifnot(identical(coef(po),numeric(0)))
coef(po,c("r","sigma")) <- 1
stopifnot(all.equal(coef(po),c(r=1,sigma=1)))
coef(po) <- NULL
coef(po,c("r","sigma"),transform=TRUE) <- 0
stopifnot(all.equal(coef(po),c(r=1,sigma=1)))
coef(po) <- NULL
coef(po) <- c(r=1,sigma=1)
stopifnot(all.equal(coef(po),c(r=1,sigma=1)))
coef(po) <- NULL
coef(po,transform=TRUE) <- c(r=0,sigma=0)
stopifnot(all.equal(coef(po),c(r=1,sigma=1)))

pompExample(ou2,envir=NULL) -> ou2
ou2[[1]] -> po
po1 <- simulate(po)

as(po,"data.frame") %>% head()
as.data.frame(po1) %>% head()

obs(po)[,1:3]
obs(po,"y2")[,1:3]
try(obs(po,c("y2","z")))

states(po)
states(po1,"x1")[,1:3]
try(states(po1,"z"))
states(po1)[,1:3]

time(po)[1:3]
time(po,t0=TRUE)[1:3]

time(po) <- 1:10
try(time(po) <- c("bob","nancy"))
time(po1,t0=TRUE) <- 0:10
try(time(po) <- 10:0)
try(time(po,t0=TRUE) <- c(4,1:10))

window(po,end=5)
window(po,start=5)
window(po,start=5,end=10)
try(window(po,start=5,end=3))
try(window(po,start=NA,end=3))
try(window(po,start=1,end=NULL))

timezero(po)
timezero(po) <- -3
try(timezero(po) <- NA)
try(timezero(po) <- c(1,2,3))
try(timezero(po) <- 20)

coef(po)
coef(po,c("alpha.3","tau"))
try(coef(po,c("alpha.3","z")))

coef(po,"alpha.3") <- 4
coef(po,"z") <- 9
coef(po)
coef(po) <- NULL
coef(po)
coef(po) <- list(a=3,b=12)

pompExample(gompertz)
gompertz -> po

coef(po)
coef(po,transform=TRUE,pars=c("r","K"))
coef(po,"sigma",transform=TRUE) <- 0
coef(po)
coef(po,c("r","K")) <- c(a=1,b=2)
coef(po,transform=TRUE) <- c(r=1,K=1)
coef(po) <- NULL
try(coef(po,transform=FALSE) <- c(5,3))
try(coef(po,transform=TRUE) <- c(5,3))
coef(po,transform=TRUE) <- c(r=1,K=1)
coef(po)
po %>%
  window(start=5,end=20) %>%
  pomp(covar=covariate_table(times=0:20,q=0:20),
    larry=3L) -> po1
as(po1,"data.frame") %>% head()

pompExample(sir,envir=NULL) %>% extract2(1) %>%
  window(end=0.12) -> po2
po2 %>% simulate(seed=4358686) %>% as.data.frame()
pomp(po2,covar=NULL)@covar
try(po2 %>% pomp(covar="bob"))
try(po2 %>% pomp(rmeasure=function(x)x))

