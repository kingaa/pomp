png(filename="bsmc2-%02d.png",res=100)
options(digits=2)

library(pomp)
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

gompertz() -> gompertz
set.seed(398585L)

time(gompertz) <- 1:10

smc <- bsmc2(gompertz,
  rprior=Csnippet("
    K = runif(0.1,1);
    r = rlnorm(log(0.2),1);
    sigma = rlnorm(log(0.1),0.5);"),
  paramnames=c("r","K","sigma"),
  Np=1000,smooth=0.05)

plot(smc,y=NA)
plot(smc,pars=c("r","K"),thin=20)
try(plot(smc,pars="bob"))
plot(smc,pars="K")
try(plot(smc,pars=NULL))
stopifnot(sum(cond.logLik(smc))==logLik(smc))
stopifnot(length(eff.sample.size(smc)) == 10)

try(bsmc2())
try(bsmc2(3L))

po <- smc
coef(po) <- NULL
try(bsmc2(smc))
try(bsmc2(po))
try(bsmc2(po,params=NULL))
try(bsmc2(po,params="yes"))
try(bsmc2(po,params=list()))
try(bsmc2(po,params=c(1,2,3)))
theta <- coef(smc)
try(bsmc2(po,params=as.list(theta)))
try(bsmc2(po,params=as.list(theta),Np=1))
try(bsmc2(smc,Np=-10))
try(bsmc2(smc,params=theta,Np=100,
  dmeasure=Csnippet("error(\"whoa nelly!\");")))
try(bsmc2(smc,Np=100,smooth=5))
try(bsmc2(smc,Np=100,smooth=NA))
try(bsmc2(smc,Np=100,smooth=-2))
try(bsmc2(smc,Np=100,smooth=NULL))
try(bsmc2(smc,Np=100,smooth=Inf))
try(bsmc2(smc,Np=100,smooth="yes"))
try(bsmc2(smc,Np=100,smooth=c(1,2)))
try(bsmc2(smc,Np=100,smooth=list(1,2)))
try(bsmc2(smc,Np=100,rprocess=NULL))
try(bsmc2(smc,Np=100,params=NULL))
try(bsmc2(smc,Np=100,dmeasure=NULL))
try(bsmc2(smc,Np=100,rprior=NULL))

theta <- coef(gompertz)
theta["K"] <- 1
try(capture.output(bsmc2(po,Np=2,params=theta,verbose=TRUE)) -> out)

smc %>% as.data.frame() %>%
  filter(.id=="posterior") %>%
  select(-.id) -> pp

gompertz %>%
  as.data.frame() %>%
  subset(select=-X) %>%
  bsmc2(
    times="time",t0=-5,
    params=coef(gompertz),
    Np=1000,smooth=0.1,
    rprior=Csnippet("
      K = runif(0.1,1);
      r = rlnorm(log(0.2),1);
      sigma = rlnorm(log(0.1),0.5);"),
    rprocess=gompertz@rprocess,
    dmeasure=gompertz@dmeasure,
    statenames=c("X"),
    paramnames=c("r","K","sigma")) -> smc4
smc4 %>% plot()

try(gompertz %>%
    as.data.frame() %>%
    subset(select=-X) %>%
    bsmc2(
      times="time",t0=-5,
      params=coef(gompertz),
      Np=1000,smooth=0.1,
      rprior=3,
      rprocess=gompertz@rprocess,
      dmeasure=gompertz@dmeasure,
      statenames=c("X"),
      paramnames=c("r","K","sigma")))

dev.off()
