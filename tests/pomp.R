options(digits=3)

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
try(data.frame(a=1:10,b=1:10) %>% pomp(covar=1:20))
try(data.frame(a=1:10,b=1:10) %>% pomp(tcovar=1:20))
try(data.frame(a=1:10,b=1:10) %>% pomp(tcovar=1:20,covar=1:20))
try(data.frame(a=1:10,b=1:10) %>%
    pomp(tcovar=1:20,covar=1:20,times=1:10))
try(data.frame(a=1:10,b=1:10) %>%
    pomp(tcovar=1:20,covar=1:20,times=1))
try(data.frame(a=1:10,b=1:10) %>%
    pomp(tcovar=1:20,covar=1:20,times=1,t0=11))
try(data.frame(a=1:10,b=1:10) %>%
    pomp(tcovar=1:20,covar=1:20,times=1,t0=NULL))
try(data.frame(a=1:10,b=1:10) %>%
    pomp(tcovar=1:20,covar=1:20,times=1,t0=NA))
try(data.frame(a=1:10,b=1:10) %>%
    pomp(tcovar=1:20,covar=1:20,times=3,t0=0))
try(data.frame(a=1:10,b=1:10) %>%
    pomp(tcovar=1:20,covar=1:20,times="c",t0=0))
try(data.frame(a=1:10,b=1:10) %>%
    pomp(tcovar=1:20,covar=1:20,times=NA,t0=0))
try(data.frame(a=1:10,b=1:10) %>%
    pomp(tcovar=1:20,covar=1:20,times=NULL,t0=0))
try(data.frame(a=1:10,b=1:10) %>%
    pomp(tcovar=1:10,covar=data.frame(c=1:10,d=1:10),times=1,t0=0))
try(data.frame(a=1:10,b=1:10) %>%
    pomp(tcovar=1,covar=data.frame(c=1:10,d=1:10),
      covarnames=c("e","f"),times=1,t0=0))
try(data.frame(a=1:10,b=1:10) %>%
    pomp(tcovar=1,covar=data.frame(c=1:10,d=1:10),
      covarnames=c("c"),times=1,t0=0))
stopifnot(data.frame(a=1:10,b=1:10) %>%
    pomp(tcovar=1,covar=data.frame(c=0:10,d=0:10),
      covarnames="d",times=1,t0=0,bob=3) %>% class() %>%
    equals("pomp"))
try(data.frame(a=1:10,b=1:10) %>%
    pomp(tcovar=1,
      covar=data.frame(c=1:10,d=1:10,d=1:10,check.names=FALSE),
      times=1,t0=0))
stopifnot(data.frame(a=1:10,b=1:10) %>%
    pomp(tcovar=1:20,covar=1:20,times="a",t0=0) %>% class() %>%
    equals("pomp"))
try(simulate())
try(NULL %>% pomp(t0=4))
try(NULL %>%
    pomp(tcovar=1:20,covar=1:20,times="a",t0=0))
try(NULL %>%
    pomp(tcovar=1:20,covar=1:20,times=1:10,t0=3))
try(NULL %>%
    pomp(tcovar=1:20,covar=1:20,times=1:10,t0=1,rinit=3))
stopifnot(NULL %>% pomp(tcovar=1:20,covar=1:20,times=1:10,t0=1)
  %>% class() %>% equals("pomp"))

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
    covar=data.frame(a=1:20,b=1:20),tcovar="a") %>% class() %>%
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
obs(window(po,start=3,end=5))
try(obs(window(po,start=3,end=5),vars="Z"))
states(po)
simulate(window(po,start=3,end=5),seed=4358686) -> po1
states(po1)
try(states(po1,vars="Z"))
time(po1)
time(po1,t0=TRUE)
try(time(po1) <- c("bob","nancy"))
time(po1) <- 4:7
obs(po1)
try(time(po1,t0=TRUE) <- c(9,4:7))
time(po1,t0=TRUE) <- c(1,4:7)
try(window(po,start="3"))
try(window(po,start=NA))
try(window(po,end=NULL))
try(window(po,end=FALSE))
try(window(po,start=5,end=2))
window(po1,start=6) %>% obs()
try(timezero(po1) <- "bob")
try(timezero(po1) <- NULL)
try(timezero(po1) <- c(3,2))
try(timezero(po1) <- list(3))
try(timezero(po1) <- -Inf)
try(timezero(po1) <- 9)
timezero(po1) <- 4

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
coef(po1) <- NULL
stopifnot(identical(coef(po1),numeric(0)))
coef(po1,c("r","sigma")) <- 1
stopifnot(all.equal(coef(po1),c(r=1,sigma=1)))
coef(po1) <- NULL
coef(po1,c("r","sigma"),transform=TRUE) <- 0
stopifnot(all.equal(coef(po1),c(r=1,sigma=1)))
coef(po1) <- NULL
coef(po1) <- c(r=1,sigma=1)
stopifnot(all.equal(coef(po1),c(r=1,sigma=1)))
coef(po1) <- NULL
coef(po1,transform=TRUE) <- c(r=0,sigma=0)
stopifnot(all.equal(coef(po1),c(r=1,sigma=1)))

pompExample(sir,envir=NULL) %>% extract2(1) %>%
  window(end=0.12) -> po2
po2 %>% simulate(seed=4358686) %>% as.data.frame()

