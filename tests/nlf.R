options(digits=3)
png(filename="nlf-%02d.png",res=100)

library(pomp2)
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
    lags=c(4,6),ti=100,tf=2000,seed=426094906L,
    rprocess=ou2@rprocess,rmeasure=ou2@rmeasure,rinit=ou2@rinit) -> m0

stopifnot(m0(0)==m0(1))

m2 <- nlf.objfun(m0,est=estnames,seed=426094906L)
plot(sapply(seq(-1,1,by=0.1),function(x)m2(c(x,0))))

library(subplex)
m2(guess) -> ll
stopifnot(logLik(m2)==-ll)
subplex(par=guess,fn=m2,control=list(reltol=1e-3)) -> out
stopifnot(out$convergence == 0)
m2(out$par)
plot(simulate(m2))
stopifnot(logLik(m2)==-out$value)

m2s <- nlf.objfun(m2,tensor=TRUE,period=10,seed=426094906L)
plot(sapply(seq(-1,1,by=0.1),function(x)m2s(c(x,0))))
subplex(par=guess,fn=m2s,control=list(reltol=1e-3)) -> out
stopifnot(out$convergence == 0)
m2s(out$par)
plot(simulate(m2s))

m2t <- nlf.objfun(m2s,tensor=FALSE,seed=426094906L,fail.value=1e9)
plot(sapply(seq(-1,1,by=0.1),function(x)m2t(c(x,0))))
subplex(par=guess,fn=m2t,control=list(reltol=1e-3)) -> out
stopifnot(out$convergence == 0)
m2t(out$par)
plot(simulate(m2t))
plot(pfilter(m2t,dmeasure=ou2@dmeasure,Np=100))
m2t(NA)

nlf.objfun(ou2,lags=c(1,2,3),ti=100,tf=500)
nlf.objfun(ou2,lags=c(1,2,3),est=estnames,ti=100,tf=500)
nlf.objfun(ou2,lags=c(1,2,3),est=NULL,ti=100,tf=500)
try(nlf.objfun(ou2,lags=c(1,2,3),est="bob",ti=100,tf=500))
try(nlf.objfun(ou2,lags=list(1,2,3),est=NULL,ti=100,tf=500))
try(nlf.objfun(ou2,lags=c(-1,2,3),est=estnames,ti=100,tf=500))
try(nlf.objfun(ou2,lags=c(),est=estnames,ti=100,tf=500))
try(nlf.objfun(ou2,lags="of course",est=estnames,ti=100,tf=500))
try(nlf.objfun(ou2,lags=NULL,est=estnames,ti=100,tf=500))
try(nlf.objfun(ou2,lags=c(1,2,3),est=estnames,ti=-10,tf=500))
try(nlf.objfun(ou2,lags=c(1,2,3),est=estnames,ti=NULL,tf=500))
try(nlf.objfun(ou2,lags=c(1,2,3),est=estnames,ti=NA,tf=500))
try(nlf.objfun(ou2,lags=c(1,2,3),est=estnames,ti=c(12,16),tf=500))
try(nlf.objfun(ou2,lags=c(1,2,3),est=estnames,ti=100,tf=NULL))
try(nlf.objfun(ou2,lags=c(1,2,3),est=estnames,ti=100,tf=NA))
try(nlf.objfun(ou2,lags=c(1,2,3),est=estnames,ti=100,tf=c(1,2,3)))
try(nlf.objfun(ou2,lags=c(1,2,3),est=estnames,ti=100,tf=list(5)))
try(nlf.objfun(ou2,lags=c(1,2,3),est=estnames,ti=100,tf=Inf))
try(nlf.objfun(ou2,lags=c(1,2,3),est=estnames,ti=100,tf="yes"))
try(nlf.objfun(ou2,lags=c(1,2,3),est=estnames,ti=100,tf=50))
try(nlf.objfun(ou2,lags=c(1,2,3),est=estnames,ti=100,tf=150))
nlf.objfun(ou2,lags=c(1,2,3),est=estnames,ti=100,tf=300)
try(nlf.objfun(ou2,lags=c(1,2,3),nrbf=2,ti=100,tf=500))
try(nlf.objfun(ou2,lags=c(1,2,3),nrbf=NA,ti=100,tf=500))
try(nlf.objfun(ou2,lags=c(1,2,3),nrbf=NULL,ti=100,tf=500))
try(nlf.objfun(ou2,lags=c(1,2,3),nrbf="4",ti=100,tf=500))
try(nlf.objfun(ou2,lags=c(1,2,3),period=Inf,ti=100,tf=500))
nlf.objfun(ou2,lags=c(1,2,3),period=-5,ti=100,tf=500)
try(nlf.objfun(ou2,lags=c(1,2,3),ti=100,tf=500,transform.data=3))
try(nlf.objfun(ou2,lags=c(1,2,3),ti=100,tf=500,transform.data=list()))
try(nlf.objfun(ou2,lags=c(1,2,3),ti=100,tf=500,transform.data="bob"))
try(nlf.objfun(ou2,lags=c(1,2,3),ti=100,tf=500,transform.data=NULL))
try(nlf.objfun(ou2,lags=c(1,2,3),ti=100,tf=500,transform.data=NA))
try(nlf.objfun(ou2,lags=c(1,2,3),ti=100,tf=500,transform.data=function(x,y)x+y))
try(nlf.objfun(ou2,lags=c(1,2,3),ti=100,tf=500,fail.value=c(10,20)))
try(nlf.objfun(ou2,lags=c(1,2,3),ti=100,tf=500,fail.value=NULL))
try(nlf.objfun(ou2,lags=c(1,2,3),ti=100,tf=500,fail.value="no"))
capture.output(nlf.objfun(ou2,lags=c(1,2,3),ti=100,tf=500,fail.value=10,verbose=TRUE)) -> out
stopifnot(sum(grepl("logql = ",out))==1)

po <- ou2
time(po) <- c(1:5,8:10)
try(po %>% nlf.objfun(lags=c(1,2,3),period=5,,ti=100,tf=500))

po <- ou2
po@data[2,15] <- NA
stopifnot(po %>% nlf.objfun(lags=c(1,2,3),ti=100,tf=500) %>% logLik() %>% is.na())

po <- ou2
stopifnot(po %>% nlf.objfun(rmeasure=NULL,lags=c(1,2,3),ti=100,tf=500) %>% logLik() %>% is.na())

try(po %>% nlf.objfun(covar=covariate_table(t=0:100,Z=100,times="t"),ti=100,tf=500))

dev.off()
