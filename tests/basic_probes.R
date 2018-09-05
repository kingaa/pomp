options(digits=3)

library(pomp)
library(magrittr)
library(reshape2)
library(dplyr)
library(tidyr)
library(ggplot2)

set.seed(255066335)

pompExample(ou2)

plist <- list(
  mean=probe.mean("y1"),
  med=probe.median("y2"),
  v=probe.var("y1"),
  sd=probe.sd("y2"),
  pd=probe.period("y1",kernel.width=5),
  probe.quantile("y2",prob=c(0.1,0.5,0.9)),
  probe.acf("y1",lags=c(1,3,6)),
  probe.acf("y2",lags=c(1,2,3),type="cor"),
  probe.ccf(c("y1","y2"),lags=c(0,1,2)),
  probe.marginal("y1",ref=obs(ou2,"y1")),
  probe.nlar("y2",lags=c(1,1,1,2,2),powers=c(1,2,3,1,2))
)

ou2 %>% probe(probes=plist,nsim=1000) -> pb
summary(pb) -> sm
stopifnot(names(sm)==c("coef","nsim","quantiles","pvals","synth.loglik"),
  logLik(pb)==sm$synth.loglik,
  length(sm$pvals)==25,length(sm$quantiles)==25)

try(probe.mean(c("y1","y2")))
try(probe.median(c("y1","y2")))
try(probe.var(c("y1","y2")))
try(probe.sd(c("y1","y2")))
try(probe.period(c("y1","y2")))
try(probe.quantile(c("y1","y2")))
try(probe.marginal(c("y1","y2")))
try(probe.ccf("y1"))
try(probe.nlar(c("y1","y2")))

try(probe.acf(c("y1","y2"),lags=c(0,1),type="cor"))
probe.acf(c("y1","y2"),lags=c(1,5),type="cor") -> f
ou2 %>% simulate() %>% obs() %>% f() -> v
stopifnot(names(v)==c("acf.1.y1", "acf.5.y1", "acf.1.y2", "acf.5.y2"))

try(ou2 %>% simulate(rmeasure=function(...) c(y=1)) %>% obs() %>% f())
try(ou2 %>% simulate(rmeasure=function(...) c(y1=NA,y2=NA)) %>% obs() %>% f())
ou2 %>% simulate(rmeasure=function(t,...) c(y1=-t,y2=t)) %>% obs() %>% f()
probe.ccf(c("y2","y1"),lags=c(0,1,2),type="cor") -> f
ou2 %>% simulate() %>% obs() %>% f()
try(ou2 %>% simulate(rmeasure=function(...) c(y=1)) %>% obs() %>% f())
try(ou2 %>% simulate(rmeasure=function(...) c(y1=NA,y2=NA)) %>% obs() %>% f())
ou2 %>% simulate(times=1:10,rmeasure=function(t,...) c(y1=-t,y2=t)) %>% obs() %>% f()

probe.marginal("y1",ref=obs(ou2,"y1"),order=6,diff=2) -> f
ou2 %>% simulate() %>% obs() %>% f()
try(ou2 %>% simulate(rmeasure=function(...) c(y=1)) %>% obs() %>% f())
ou2 %>% simulate(rmeasure=function(...) c(y1=NA,y2=NA)) %>% obs() %>% f() -> x
stopifnot(all(is.na(x)))
try(ou2 %>% simulate(times=1:10,rmeasure=function(t,...) c(y1=-t,y2=t)) %>% obs() %>% f())
ou2 %>% simulate(rmeasure=function(t,...) c(y1=-t,y2=t)) %>% obs() %>% f() -> x
stopifnot(x==0)

try(probe.nlar("y1",lags=c(0,-1)))
try(probe.nlar("y1",lags=c(0,-1),powers=2))
try(probe.nlar("y1",lags=c(0,1,2),powers=c(0,-1)))
try(probe.nlar("y1",lags=c(0,1,NA),powers=c(0,1)))
try(probe.nlar("y1",lags=c(0,1,2),powers=NA))
try(probe.nlar("y1",lags=c(0,1,2),powers=NULL))
try(probe.nlar("y1",lags=list(0,1,2),powers=1))
try(probe.nlar("y1",lags=list(1,2),powers=c(1,2,3)))
try(probe.nlar("y1",lags=list(1,2,3),powers=c(1,2)))

plist <- list(
  probe.nlar("y1",lags=c(1,2,3),powers=2),
  probe.nlar("y2",lags=1,powers=c(1,2,3))
)
probe.nlar("y2",lags=1,powers=c(1,2,3)) -> f
ou2 %>% simulate() %>% obs() %>% f()
try(ou2 %>% simulate(rmeasure=function(...) c(y=1)) %>% obs() %>% f())
ou2 %>% simulate(rmeasure=function(...) c(y1=NA,y2=NA)) %>% obs() %>% f() -> x
stopifnot(x==0)
ou2 %>% simulate(times=1:10,rmeasure=function(t,...) c(y1=-t,y2=t)) %>% obs() %>% f()
