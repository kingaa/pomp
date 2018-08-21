options(digits=3)
png(filename="spect_match-%02d.png",res=100)

library(pomp)
library(magrittr)

pompExample(ou2)

ou2 %>%
  as.data.frame() %>%
  spect.match.objfun(
    times="time",t0=0,
    rinit=ou2@rinit,
    rprocess=ou2@rprocess,
    rmeasure=ou2@rmeasure,
    kernel.width=3,
    params=coef(ou2),
    nsim=100,
    seed=5069977
    ) -> f

f(0)

f %>% spect.match.objfun(est=c("alpha.1")) -> f1
plot(sapply(seq(0.3,1.2,by=0.1),f1),log='y')

f1(1.1)
plot(spect(f1))
library(subplex)
subplex(fn=f1,par=0.4,control=list(reltol=1e-3)) -> out
f1(out$par)
plot(spect(f1))

dev.off()
