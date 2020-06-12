options(digits=3)
png(filename="wpfilter-%02d.png",res=100)

library(pomp)

set.seed(9994847L)

ricker() -> po
try(wpfilter(po))
try(wpfilter(po,Np=NULL))
try(wpfilter(po,Np=-10))
try(wpfilter(po,Np=c(10,20,30)))

wpfilter(po,Np=200,trigger=1,target=1) %>% plot()
wpfilter(po,Np=200,trigger=1,target=0.5) %>% plot()
wpfilter(po,Np=200,trigger=1,target=0) %>% plot()
wpfilter(po,Np=200,trigger=0.8) %>% plot()

wpfilter(po,Np=1000,trigger=0.5,target=0.5) %>% logLik()
wpfilter(po,Np=1000,trigger=0.5,target=0.2) %>% logLik()
wpfilter(po,Np=1000,trigger=0.5,target=0) %>% logLik()
wpfilter(po,Np=1000,trigger=1,target=0) %>% logLik()
pfilter(po,Np=1000) %>% logLik()

do.call(c,replicate(n=10,pfilter(window(po,end=4),Np=1000))) %>% logLik() %>% logmeanexp(se=TRUE)
do.call(c,replicate(n=10,wpfilter(window(po,end=4),Np=1000))) %>% logLik() %>% logmeanexp(se=TRUE)

do.call(c,replicate(n=10,pfilter(window(po,end=20),Np=1000))) %>% logLik() %>% logmeanexp(se=TRUE)
do.call(c,replicate(n=10,wpfilter(window(po,end=20),Np=1000,trigger=0.8))) %>% logLik() %>% logmeanexp(se=TRUE)

set.seed(9994847L)
try(wpfilter())
try(wpfilter("bob"))
try(wpfilter(list(3,2,1)))
wpfilter(po,Np=100,dmeasure=function(...,log)-Inf)

wpfilter(po,Np=100) %>% wpfilter() -> pf
wpfilter(po,Np=100) %>% wpfilter(target=0,trigger=0.1,Np=200) -> pf
try(wpfilter(po,dmeasure=NULL))
try(wpfilter(po,rprocess=NULL))
try(wpfilter(po,Np=100,trigger=-1))
try(wpfilter(po,Np=100,trigger=NULL))
try(wpfilter(po,Np=100,trigger=c(0,1)))
try(wpfilter(po,Np=100,trigger=NA))
try(wpfilter(po,Np=100,target=-1))
try(wpfilter(po,Np=100,target=2))
try(wpfilter(po,Np=100,target=NULL))
try(wpfilter(po,Np=100,target=c(0,1)))
try(wpfilter(po,Np=100,target=NaN))

try(wpfilter(po,Np=100,dmeasure=function(...,log)sample(c(0,Inf),1)))

po %>%
  as.data.frame() %>%
  wpfilter(
    times="time",t0=0,Np=100,
    params=c(N0=7,phi=10,r=40),
    rinit=function(N0,...)c(N=N0),
    dmeasure=function(y,N,phi,...,log) dnorm(x=y,mean=phi*N,sd=sqrt(phi*N),log=log),
    rprocess=discrete_time(function(r,N,...) c(N=r*N*exp(-N+rnorm(1))),delta.t=1)
  ) %>% plot()

dev.off()
