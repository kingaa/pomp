options(digits=3)
png(filename="traj_match-%02d.png",res=100)

library(pomp2)
library(magrittr)

pompExample(ou2)

try(traj.match.objfun())
try(traj.match.objfun("bob"))
try(ou2 %>% as.data.frame() %>% traj.match.objfun())
try(ou2 %>% as.data.frame() %>% traj.match.objfun(times="time",t0=0))

ou2 %>%
  as.data.frame() %>%
  traj.match.objfun(
    times="time",t0=0,
    rinit=ou2@rinit,
    skeleton=ou2@skeleton,
    dmeasure=ou2@dmeasure,
    params=coef(ou2)
  ) -> f

stopifnot(f(0)==f(1))

f %>% traj.match.objfun(est=c("alpha.1")) -> f1
plot(sapply(seq(0.1,0.9,by=0.1),f1),xlab="",ylab="")

f1(1.1)
matplot(t(trajectory(f1)[,1,]),type="l",ylab="y")
library(subplex)
subplex(fn=f1,par=0.4,control=list(reltol=1e-3)) -> out
f1(out$par)

try(traj.match.objfun(f1,est="harry"))

f1 %>% as("pomp")
f1 %>% as("data.frame") %>% names()

f1 %>% traj.match.objfun(fail.value=1e10) -> f2
f2(NA)

dev.off()
