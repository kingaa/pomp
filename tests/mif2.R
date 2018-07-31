options(digits=3)
png(filename="mif2-%02d.png",res=100)

set.seed(857075216L)

library(pomp)
library(dplyr)
library(magrittr)

pompExample(gompertz,envir=NULL) -> po
po[[1]] %>% window(end=10) -> po

mif2(po,Nmif=50,Np=100,transform=TRUE,cooling.fraction.50=0.5,
  rw.sd=rw.sd(sigma=0.02,K=0.02,r=0.02)) -> mf
plot(mf)

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
try(mif2(po,Nmif=1,Np=c(rep(100,11),40),rw.sd=rw.sd(sigma=0.1),cooling.frac=0.5))
try(mif2(po,Nmif=1,Np=100,rw.sd=3))
try(mif2(po,Nmif=1,Np=100,rw.sd=rw.sd()))
try(mif2(po,Nmif=1,Np=100,rw.sd=rw.sd(a=9)))
try(mif2(po,Nmif=1,Np=100,rw.sd=rw.sd(sigma=0.1)))
try(mif2(po,Nmif=1,Np=100,rw.sd=rw.sd(NULL)))
try(mif2(po,Nmif=1,Np=100,rw.sd=rw.sd(sigma=0.1),cooling.fraction.50=12))
try(mif2(po,Nmif=1,Np=100,rw.sd=rw.sd(sigma=0.1),cooling.fraction.50=NA))
try(mif2(po,Nmif=1,Np=100,rw.sd=rw.sd(sigma=0.1),cooling.fraction.50=c(3,2)))
try(mif2(po,Nmif=1,Np=100,rw.sd=rw.sd(sigma=0.1),cooling.fraction.50=NULL))
try(mif2(po,Nmif=1,Np=100,rw.sd=rw.sd(sigma=0.1),cooling.fraction.50=0.1,
  cooling.type="quick"))
try(mif2(po,Nmif=1,Np=100,rw.sd=rw.sd(sigma=0.1),cooling.fraction.50=0.1,
  cooling.type="geometric",tol=-3))
try(mif2(po,start=NULL,Nmif=1,Np=100))
try(mif2(po,start=list(),Nmif=1,Np=100))
try(mif2(po,start=list(NULL),Nmif=1,Np=100))
try(mif2(po,start=c(3,2,1),Nmif=1,Np=100))
try(mif2(po,Nmif=1,Np=100:1000,rw.sd=rw.sd(sigma=0.1)))
mif2(po,Nmif=3,Np=100,rw.sd=rw.sd(sigma=0.01,X.0=ivp(0.01)),
  cooling.fraction.50=0.1,cooling.type="geometric",tol=1e-10,
  transform=TRUE,start=as.list(coef(po)))

theta <- coef(po)
theta["sigma"] <- 0.2
po %>%
  pfilter(Np=100,params=theta) %>%
  mif2(Nmif=3,rw.sd=rw.sd(sigma=0.01,X.0=ivp(0.01)),
    cooling.fraction.50=0.5,transform=TRUE) %>%
  mif2() %>% continue(Nmif=3,cooling.fraction.50=0.1) %>% plot()

dev.off()
