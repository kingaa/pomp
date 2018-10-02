options(digits=3)
png(filename="probe_match-%02d.png",res=100)

library(pomp2)

gompertz() -> po
po %>% as.data.frame() %>% subset(select=-X) -> dat

plist <- list(
  mean=probe.mean("Y",trim=0.1,transform=sqrt),
  sd=probe.sd("Y",transform=sqrt),
  probe.marginal("Y",ref=obs(po)),
  probe.acf("Y",lags=c(1,3,5),type="correlation",transform=sqrt),
  probe.quantile("Y",prob=c(0.25,0.75),na.rm=TRUE)
)

try(dat %>% probe_objfun())
try(dat %>% probe_objfun(times="time",t0=0))

dat %>%
  probe_objfun(
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

f %>% probe_objfun(est=c("K"),seed=580656309) -> f1
plot(sapply(seq(0.8,1.6,by=0.1),f1))

f1(1.1)
library(subplex)
subplex(fn=f1,par=0.4,control=list(reltol=1e-3)) -> out
f1(out$par)

try(probe_objfun())
try(probe_objfun("bob"))

try(probe_objfun(f,est="harry"))

f1 %>% as("probed_pomp") %>% plot()

f1 %>% summary() %>% names()

f1 %>% probe() %>% plot()

f1 %>% as.pomp() %>% as.data.frame() %>% names()

f1 %>% probe_objfun(fail.value=1e10) -> f2

f1 %>% spect(kernel.width=3,nsim=100,seed=748682047) %>% plot()

f1 %>% as("pomp")
f1 %>% as("data.frame") %>% names()

po %>% probe_objfun(nsim=100,probes=function(x)1,fail.value=1e9) -> f2
logLik(f2)
f2(1)

dev.off()
