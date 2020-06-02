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
po %>% wpfilter(Np=100) -> pf
plot(pf)

dev.off()
