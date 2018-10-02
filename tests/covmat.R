options(digits=3)

library(pomp2)
library(tidyr)
library(dplyr)

gompertz() -> gompertz

set.seed(1546383977L)

plist <- list(
  y1.mean=probe.mean(var="Y"),
  probe.acf(var="Y",lags=c(0,4,8))
)

gompertz %>%
  probe(probes=plist,nsim=100) %>%
  as.data.frame() %>%
  filter(.id=="sim") %>%
  select(-.id) %>%
  gather(variable,value) %>%
  group_by(variable) %>%
  summarize(value=(sd(value))) %>%
  ungroup() %>%
  spread(variable,value) %>%
  unlist() -> scale.dat

abc(
  gompertz,
  dprior=Csnippet("
    lik = dunif(K,0,2,1)+dunif(r,0,1,1)+
      dunif(tau,0,1,1)+dunif(sigma,0,1,1);
    lik = (give_log) ? lik : exp(lik);"),
  paramnames=c("K","sigma","tau","r"),
  Nabc=100,probes=plist,
  scale=scale.dat,epsilon=10,
  proposal=mvn.diag.rw(c(K=0.01,r=0.01,sigma=0.01,tau=0.01))
) -> a1

replicate(3,
  abc(
    a1,
    Nabc=500,probes=plist,
    scale=scale.dat,epsilon=10,
    proposal=mvn.diag.rw(c(K=0.01,r=0.01,sigma=0.01,tau=0.01))
  )) -> a1
do.call(c,a1) -> a1

covmat(a1[[1]]) -> v1
covmat(a1) -> v2
covmat(a1,thin=20) -> v3
stopifnot(
  dim(v1)==dim(v2),
  dim(v1)==dim(v3),
  identical(dimnames(v1),dimnames(v2)),
  identical(dimnames(v1),dimnames(v3))
)

po <- window(gompertz,end=10)

prop1 <- mvn.diag.rw(c(r=0.01,sigma=0.01))

mcmc1 <- pmcmc(po,Nmcmc=100,Np=100,dprior=Csnippet("
    lik = dunif(r,0,1,1)+dnorm(sigma,0,1,1);
    lik = (give_log) ? lik : exp(lik);"),
  paramnames=c("r","sigma"),
  proposal=prop1)

covmat(mcmc1) -> v1
covmat(c(mcmc1,mcmc1)) -> v2
stopifnot(
  dim(v1)==dim(v2),
  identical(dimnames(v1),dimnames(v2))
)
