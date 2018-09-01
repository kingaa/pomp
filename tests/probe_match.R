options(digits=3)
png(filename="probe_match-%02d.png",res=100)

library(pomp)
library(magrittr)

pompExample(gompertz)
gompertz -> po

plist <- list(
  mean=probe.mean("Y",trim=0.1,transform=sqrt),
  sd=probe.sd("Y",transform=sqrt),
  probe.marginal("Y",ref=obs(po)),
  probe.acf("Y",lags=c(1,3,5),type="correlation",transform=sqrt),
  probe.quantile("Y",prob=c(0.25,0.75),na.rm=TRUE)
)

try(gompertz %>% as.data.frame() %>% probe.match.objfun())

gompertz %>%
  as.data.frame() %>%
  probe.match.objfun(
    times="time",t0=0,
    rinit=po@rinit,
    rprocess=po@rprocess,
    rmeasure=po@rmeasure,
    probes=plist,
    params=coef(po),
    nsim=100,
    seed=5069977
    ) -> f

stopifnot(f(0)==f(1))

f %>% probe.match.objfun(est=c("K"),seed=580656309) -> f1
plot(sapply(seq(0.8,1.6,by=0.1),f1))

f1(1.1)
library(subplex)
subplex(fn=f1,par=0.4,control=list(reltol=1e-3)) -> out
f1(out$par)

try(probe.match.objfun())
try(probe.match.objfun("bob"))

try(probe.match.objfun(f,est="harry"))

f1 %>% as("probed_pomp") %>% plot()

f1 %>% summary() %>% names()

f1 %>% probe() %>% plot()

f1 %>% as.pomp() %>% as.data.frame() %>% names()

f1 %>% probe.match.objfun(fail.value=1e10) -> f2

dev.off()
