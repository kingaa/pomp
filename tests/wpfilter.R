options(digits=3)
png(filename="wpfilter-%02d.png",res=100)

library(pomp)

set.seed(9994847L)

ricker() -> po
try(wpfilter(po))
try(wpfilter(po,Np=NULL))
try(wpfilter(po,Np=-10))
try(wpfilter(po,Np=c(10,20,30)))
wpfilter(po,Np=ceiling(runif(52,min=10,max=100)))

wpfilter(po,Np=200,trigger=1,target=1) %>% plot()
wpfilter(po,Np=200,trigger=1,target=0.5) %>% plot()
wpfilter(po,Np=200,trigger=1,target=0) %>% plot()
wpfilter(po,Np=200,trigger=0.8) %>% plot()

wpfilter(po,Np=1000,trigger=0.5,target=0.5) %>% logLik()
wpfilter(po,Np=1000,trigger=0.5,target=0.2) %>% logLik()
wpfilter(po,Np=1000,trigger=0.5,target=0) %>% logLik()
wpfilter(po,Np=1000,trigger=1,target=0) %>% logLik()

do.call(c,replicate(n=10,pfilter(window(po,end=4),Np=1000))) %>% logLik() %>% logmeanexp(se=TRUE)
do.call(c,replicate(n=10,wpfilter(window(po,end=4),Np=1000))) %>% logLik() %>% logmeanexp(se=TRUE)

do.call(c,replicate(n=10,pfilter(window(po,end=20),Np=1000))) %>% logLik() %>% logmeanexp(se=TRUE)
do.call(c,replicate(n=10,wpfilter(window(po,end=20),Np=1000,trigger=0.8))) %>% logLik() %>% logmeanexp(se=TRUE)

dev.off()
