options(digits=3)
png(filename="traj_match-%02d.png",res=100)

library(pomp)

ou2() -> ou2

try(traj_objfun())
try(traj_objfun("bob"))
try(ou2 %>% as.data.frame() %>% traj_objfun())
try(ou2 %>% as.data.frame() %>% traj_objfun(times="time",t0=0))

ou2 %>%
  as.data.frame() %>%
  traj_objfun(
    times="time",t0=0,
    rinit=ou2@rinit,
    skeleton=ou2@skeleton,
    dmeasure=ou2@dmeasure,
    params=coef(ou2)
  ) -> f

stopifnot(f(0)==f(1))

f %>% traj_objfun(est=c("alpha_1")) -> f1
plot(sapply(seq(0.1,0.9,by=0.1),f1),xlab="",ylab="")

f1(1.1)
matplot(t(trajectory(f1)[,1,]),type="l",ylab="y")
library(subplex)
subplex(fn=f1,par=0.4,control=list(reltol=1e-3)) -> out
f1(out$par)

try(traj_objfun(f1,est="harry"))

f1 %>% as("pomp")
f1 %>% as("data.frame") %>% names()

f1 %>% traj_objfun(fail.value=1e10) -> f2
f2(NA)

## ------------------
## cf issue #149
## ------------------

enames <- c("gamma","iota","rho","k")

sir() |>
  traj_objfun(
    dmeasure=Csnippet("
      lik = dnbinom_mu(nearbyint(reports),1/k,rho*cases,give_log);
    "),
    rmeasure=NULL,
    params=c(coef(sir()),k=1),
    partrans=parameter_trans(
      log=c("gamma","iota","k"),
      logit="rho"
    ),
    est=enames,
    paramnames=enames,
    statenames="cases"
  ) -> ofun

theta <- c(gamma=10,iota=1,rho=0.2,k=1)

subplex(
  fn=ofun,
  par=partrans(ofun,theta,dir="toEst")
) -> fit

invisible(ofun(fit$par))

library(ggplot2)

ofun %>%
  trajectory(format="d") %>%
  ggplot(aes(x=time,y=coef(ofun,"rho")*cases))+
  geom_line()+
  geom_point(
    data=as(ofun,"data.frame"),
    mapping=aes(x=time,y=reports)
  )+
  theme_bw()
## ------------------

dev.off()
