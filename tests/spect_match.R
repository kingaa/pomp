options(digits=3)
png(filename="spect_match-%02d.png",res=100)

library(pomp2)
library(magrittr)

ou2() -> ou2
ou2 %>% as.data.frame() %>% subset(select=c(time,y1,y2)) -> dat

try(dat %>% spect.objfun())
try(dat %>% spect.objfun(times="time",t0=0))

dat %>%
  spect.objfun(
    times="time",t0=0,
    rinit=ou2@rinit,
    rprocess=ou2@rprocess,
    rmeasure=ou2@rmeasure,
    kernel.width=3,
    params=coef(ou2),
    nsim=100,
    seed=5069977
    ) -> f

stopifnot(f(0)==f(1))
stopifnot(logLik(f)==-f(0))

f %>% spect.objfun(est=c("alpha_1"),seed=580656309) -> f1
plot(sapply(seq(0.3,1.2,by=0.1),f1),log='y')

f1(1.1)
plot(spect(f1))
library(subplex)
subplex(fn=f1,par=0.4,control=list(reltol=1e-3)) -> out
f1(out$par)

try(spect.objfun())
try(spect.objfun("bob"))
try(spect.objfun(f1,est="harry"))

f1 %>% as("spectd_pomp") %>% plot()

f1 %>% summary() %>% names()

f1 %>% spect() %>% plot()

f1 %>% as("pomp")
f1 %>% as("data.frame") %>% names()

f1 %>% spect.objfun(fail.value=1e10) -> f2

try(spect.objfun(f2,weights="heavy"))
try(spect.objfun(f2,weights=c(3,4,5)))
spect.objfun(f2,weights=exp(-seq(0,1,length=50)))
try(spect.objfun(f2,weights=function(f)1-4*f))
try(spect.objfun(f2,weights=function(f)stop("oh no!")))
spect.objfun(f2,seed=5069977,weights=function(f)exp(-f/0.1)) -> f2
subplex(fn=f2,par=out$par,control=list(reltol=1e-3)) -> out
f2(out$par)
summary(f2) %>% names()
plot(spect(f2))

f2 %>% probe(nsim=100,seed=639095851,
  probes=list(probe.mean("y1"),probe.mean("y2"))) %>%
  plot()

dev.off()
