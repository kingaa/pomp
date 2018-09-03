options(digits=3)
png(filename="nlf-%02d.png",res=100)

library(pomp)
library(magrittr)

set.seed(583615606L)

pompExample(ou2)
estnames=c("alpha.2","alpha.3")
ou2 %>% window(end=40) %>% as.data.frame() -> dat
guess <- coef(ou2,estnames,transform=TRUE)

try(nlf.objfun())
try(nlf.objfun("bob"))

dat %>%
  nlf.objfun(times="time",t0=0,params=coef(ou2),
    lags=c(4,6),ntransient=100,nasymp=2000,seed=426094906L,
    rprocess=ou2@rprocess,rmeasure=ou2@rmeasure,rinit=ou2@rinit) -> m0

stopifnot(m0(0)==m0(1))

m2 <- nlf.objfun(m0,est=estnames,seed=426094906L)

library(subplex)
subplex(par=guess,fn=m2) -> out
stopifnot(out$convergence == 0)
m2(out$par)
plot(simulate(m2))

m2s <- nlf.objfun(m2,tensor=TRUE,period=10,seed=426094906L)
subplex(par=guess,fn=m2s) -> out
stopifnot(out$convergence == 0)
m2s(out$par)
plot(simulate(m2s))

m2t <- nlf.objfun(m2s,tensor=FALSE,seed=426094906L)
subplex(par=guess,fn=m2t) -> out
stopifnot(out$convergence == 0)
m2t(out$par)
plot(simulate(m2t))

nlf.objfun(ou2,lags=c(1,2,3),nasymp=100,ntransient=100)
nlf.objfun(ou2,lags=c(1,2,3),est=estnames,nasymp=100,ntransient=100)
nlf.objfun(ou2,lags=c(1,2,3),est=NULL,nasymp=100,ntransient=100)
try(nlf.objfun(ou2,lags=list(1,2,3),est=NULL,nasymp=100,ntransient=100))
try(nlf.objfun(ou2,lags=c(-1,2,3),est=estnames,nasymp=100,ntransient=100))
try(nlf.objfun(ou2,lags=c(),est=estnames,nasymp=100,ntransient=100))
try(nlf.objfun(ou2,lags="of course",est=estnames,nasymp=100,ntransient=100))
try(nlf.objfun(ou2,lags=NULL,est=estnames,nasymp=100,ntransient=100))
try(nlf.objfun(ou2,lags=c(1,2,3),est=estnames,nasymp=-100,ntransient=100))
try(nlf.objfun(ou2,lags=c(1,2,3),est=estnames,nasymp=NULL,ntransient=100))
try(nlf.objfun(ou2,lags=c(1,2,3),est=estnames,nasymp=NA,ntransient=100))
try(nlf.objfun(ou2,lags=c(1,2,3),est=estnames,nasymp=c(12,16),ntransient=100))
try(nlf.objfun(ou2,lags=c(1,2,3),est=estnames,nasymp=100,ntransient=NULL))
try(nlf.objfun(ou2,lags=c(1,2,3),est=estnames,nasymp=100,ntransient=NA))
try(nlf.objfun(ou2,lags=c(1,2,3),est=estnames,nasymp=100,ntransient=c(1,2,3)))
try(nlf.objfun(ou2,lags=c(1,2,3),est=estnames,nasymp=100,ntransient=list(5)))
try(nlf.objfun(ou2,lags=c(1,2,3),est=estnames,nasymp=100,ntransient=Inf))
try(nlf.objfun(ou2,lags=c(1,2,3),est=estnames,nasymp=100,ntransient="yes"))
try(nlf.objfun(ou2,lags=c(1,2,3),est=estnames,nasymp=100,ntransient=3))
try(nlf.objfun(ou2,lags=c(1,2,3),nrbf=2,nasymp=100,ntransient=300))
try(nlf.objfun(ou2,lags=c(1,2,3),nrbf=NA,nasymp=100,ntransient=300))
try(nlf.objfun(ou2,lags=c(1,2,3),nrbf=NULL,nasymp=100,ntransient=300))
try(nlf.objfun(ou2,lags=c(1,2,3),nrbf="4",nasymp=100,ntransient=300))
try(nlf.objfun(ou2,lags=c(1,2,3),period=Inf,nasymp=100,ntransient=300))
nlf.objfun(ou2,lags=c(1,2,3),period=-5,nasymp=100,ntransient=300)
try(nlf.objfun(ou2,lags=c(1,2,3),nasymp=100,ntransient=300,transform.data=3))
try(nlf.objfun(ou2,lags=c(1,2,3),nasymp=100,ntransient=300,transform.data=list()))
try(nlf.objfun(ou2,lags=c(1,2,3),nasymp=100,ntransient=300,transform.data="bob"))
try(nlf.objfun(ou2,lags=c(1,2,3),nasymp=100,ntransient=300,transform.data=NULL))
try(nlf.objfun(ou2,lags=c(1,2,3),nasymp=100,ntransient=300,transform.data=NA))
try(nlf.objfun(ou2,lags=c(1,2,3),nasymp=100,ntransient=300,fail.value=c(10,20)))
try(nlf.objfun(ou2,lags=c(1,2,3),nasymp=100,ntransient=300,fail.value=NULL))
try(nlf.objfun(ou2,lags=c(1,2,3),nasymp=100,ntransient=300,fail.value="no"))
nlf.objfun(ou2,lags=c(1,2,3),nasymp=100,ntransient=300,fail.value=10)

po <- ou2
time(po) <- c(1:5,8:10)
try(nlf.objfun(po,lags=c(1,2,3),period=-5,nasymp=100,ntransient=300))

dev.off()
