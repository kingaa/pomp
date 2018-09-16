options(digits=3)
png(filename="mif2-%02d.png",res=100)

set.seed(857075216L)

library(pomp2)
library(dplyr)
library(reshape2)
library(magrittr)

gompertz() %>% window(end=10) -> po

mif2(po,Nmif=50,Np=100,cooling.fraction.50=0.5,
  rw.sd=rw.sd(sigma=0.02,K=0.02,r=0.02)) -> mf1
mif2(po,Nmif=50,Np=100,cooling.fraction.50=0.5,
  rw.sd=rw.sd(sigma=0.02,K=0.02,r=0.02)) -> mf2
plot(mf1)
plot(c(a=mf1,b=mf2) -> mfl,y=NA)
c(a=mf1,b=c(mf1,mf2))
mfl[1]
mfl["b"]
mfl[5]
traces(mfl) -> tr
stopifnot(
  length(tr)==2,
  names(tr)==c("a","b"),
  dim(tr$a)==dim(tr$b),
  identical(dimnames(tr$a),dimnames(tr$b)),
  colnames(tr$a)==c("loglik","nfail","K","r","sigma","tau","X_0")
)
coef(mfl) -> m
stopifnot(
  dim(m)==c(5,2),
  colnames(m)==c("a","b"),
  rownames(m)==names(coef(po))
)

mfl %>%
  traces(transform=TRUE,pars=c("r","sigma")) -> tr
stopifnot(
  length(tr)==2,
  names(tr)==c("a","b"),
  dim(tr$a)==dim(tr$b),
  identical(dimnames(tr$a),dimnames(tr$b)),
  colnames(tr$a)==c("r","sigma")
)


try(mfl %>% traces(pars="bob"))

try(mif2())
try(mif2("po"))
try(mif2(po,Nmif=NA,Np=100))
try(mif2(po,Nmif=NULL,Np=100))
try(mif2(po,Nmif=-10,Np=100))
try(mif2(po,Nmif=c(10,20),Np=100))
try(mif2(po,Nmif=1,Np=function(k)c(10,20)))
try(mif2(po,Nmif=1,Np="bob"))
try(mif2(po,Nmif=list(),Np=100))
try(mif2(po,Nmif=1,Np=Inf))
try(mif2(po,Nmif=1,Np=100))
try(mif2(po,Nmif=1,Np=NULL))
try(mif2(po,Nmif=1,Np=c(3,4)))
mif2(po,Nmif=1,Np=c(rep(100,11),40),rw.sd=rw.sd(sigma=0.1),cooling.frac=0.5)
try(mif2(po,Nmif=1,Np=100,rw.sd=3))
try(mif2(po,Nmif=1,Np=100,rw.sd=rw.sd()))
try(mif2(po,Nmif=1,Np=100,rw.sd=rw.sd(a=9)))
try(mif2(po,Nmif=1,Np=100,rw.sd=rw.sd(sigma=1:1000)))
try(mif2(po,Nmif=1,Np=100,rw.sd=rw.sd(sigma=0.1)))
try(mif2(po,Nmif=1,Np=100,rw.sd=rw.sd(NULL)))
try(mif2(po,Nmif=1,Np=100,rw.sd=rw.sd(sigma=0.1),cooling.fraction.50=12))
try(mif2(po,Nmif=1,Np=100,rw.sd=rw.sd(sigma=0.1),cooling.fraction.50=NA))
try(mif2(po,Nmif=1,Np=100,rw.sd=rw.sd(sigma=0.1),cooling.fraction.50=c(0.1,1)))
try(mif2(po,Nmif=1,Np=100,rw.sd=rw.sd(sigma=0.1),cooling.fraction.50=NULL))
try(mif2(po,Nmif=1,Np=100,rw.sd=rw.sd(sigma=0.1),cooling.fraction.50=0.1,
  cooling.type="geometric",tol=-3))
try(mif2(po,params=NULL,Nmif=1,Np=100,rw.sd=rw.sd(sigma=0.1)))
try(mif2(po,params=list(),Nmif=1,Np=100,rw.sd=rw.sd(sigma=0.1)))
try(mif2(po,params=list(NULL),Nmif=1,Np=100,rw.sd=rw.sd(sigma=0.1)))
try(mif2(po,params=c(3,2,1),Nmif=1,Np=100,rw.sd=rw.sd(sigma=0.1)))
try(mif2(po,Nmif=1,Np=100:1000,rw.sd=rw.sd(sigma=0.1)))
mif2(po,Nmif=2,Np=50,rw.sd=rw.sd(sigma=0.01,X_0=ivp(0.01)),
  cooling.fraction.50=0.1,cooling.type="geometric",tol=1e-10,
  params=as.list(coef(po)))
try(mif2(po,Nmif=2,Np=100,rw.sd=rw.sd(sigma=0.01,X_0=ivp(0.01)),
  cooling.fraction.50=0.1,rprocess=onestep(function(x,t,params,covars,delta.t,...)stop("boink"))))
try(mif2(po,Nmif=2,Np=100,rw.sd=rw.sd(sigma=0.01,X_0=ivp(0.01)),
  cooling.fraction.50=0.1,dmeasure=function(log,...)stop("blop")))
try(mif2(po,Nmif=2,Np=100,rw.sd=rw.sd(sigma=0.01,X_0=ivp(0.01)),
  cooling.fraction.50=0.1,dmeasure=function(log,...)NA))
mif2(po,Nmif=2,Np=50,rw.sd=rw.sd(sigma=0.01),cooling.fraction.50=0.1,
  drpocess="oops",
  dmeasure=function(log,...)0) -> mf3
try(mif2(mf3,max.fail=1))
try(mif2(po,Nmif=2,Np=50,rw.sd=rw.sd(sigma=0.01),cooling.fraction.50=0.1,dmeasure=NULL))
try(mif2(po,Nmif=2,Np=50,rw.sd=rw.sd(sigma=0.01),cooling.fraction.50=0.1,rprocess=NULL))
try(mif2(po,Nmif=2,Np=50,rw.sd=rw.sd(sigma=0.01),cooling.fraction.50=0.1,params=NULL))

theta <- coef(po)
theta["sigma"] <- 0.2
po %>%
  pfilter(Np=100,params=theta) %>%
  mif2(Nmif=3,rw.sd=rw.sd(sigma=0.01,X_0=ivp(0.01)),
    cooling.fraction.50=0.5) %>%
  mif2() %>% continue(Nmif=3,cooling.fraction.50=0.1) %>% plot()

capture.output(
  mif2(po,Nmif=2,Np=100,rw.sd=rw.sd(sigma=0.01,X_0=ivp(0.01)),
    cooling.fraction.50=1,cooling.type="hyperbolic",tol=1e-10,
    params=as.list(coef(po)),verbose=TRUE),
  type="output"
) -> out
stopifnot(sum(grepl("mif2 pfilter timestep",out))==4,
  sum(grepl("mif2 iteration",out))==2)
capture.output(
  mif2(po,Nmif=2,Np=100,rw.sd=rw.sd(sigma=0.01,X_0=ivp(0.01)),
    cooling.fraction.50=1,cooling.type="hyperbolic",tol=10,
    params=as.list(coef(po)),verbose=TRUE),
  type="message"
) -> out
stopifnot(sum(grepl("filtering failure at time",out))==20)

po %>%
  as.data.frame() %>%
  subset(select=-X) %>%
  mif2(Nmif=3,Np=100,
    times="time",t0=0,
    params=c(sigma=5),
    rw.sd=rw.sd(sigma=0.01),
    cooling.fraction.50=1,cooling.type="hyperbolic",
    rprocess=onestep(function(X,...)c(X=X)),
    dmeasure=function(Y,X,sigma,log,...)dnorm(x=Y,mean=X,sd=sigma,log=log),
    rinit=function(...)c(X=0)
  )

dev.off()
