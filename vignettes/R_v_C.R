## ----packages-----------------------------------------------------------------
library(pomp)
library(ggplot2)
library(tidyr)


## ----seed,echo=FALSE----------------------------------------------------------
set.seed(56300069)


## ----R1-----------------------------------------------------------------------
simulate(times=1:100,t0=0,
  params=c(K=1,r=0.1,sigma=0.1,tau=0.1,X.0=1),
  rprocess=discrete_time(
    step.fun=function (X,r,K,sigma,...,delta.t) {
      eps <- exp(rnorm(n=1,mean=0,sd=sigma))
      S <- exp(-r*delta.t)
      c(X=K^(1-S)*X^S*eps)
    },
    delta.t=1 
  ),
  rmeasure=function (X, tau, ...) {
    c(Y=rlnorm(n=1,meanlog=log(X),sdlog=tau))
  },
  dmeasure=function (tau, X, Y, ..., log) {
    dlnorm(x=Y,meanlog=log(X),sdlog=tau,log=log)
  }
) -> gompertz


## ----R2-----------------------------------------------------------------------
gompertz |>
  as.data.frame() |>
  pivot_longer(c(X,Y)) |>
  ggplot(aes(x=time,y=value,color=name))+
  geom_line()+
  labs(y="X, Y")+
  theme_bw()


## ----C1-----------------------------------------------------------------------
simulate(
  times=0:100,t0=0,
  params=c(K=1,r=0.1,sigma=0.1,tau=0.1,X.0=1),
  dmeasure=Csnippet("
    lik = dlnorm(Y,log(X),tau,give_log);"
  ),
  rmeasure=Csnippet("
    Y = rlnorm(log(X),tau);"
  ),
  rprocess=discrete_time(
    step.fun=Csnippet("
    double S = exp(-r*dt);
    double logeps = (sigma > 0.0) ? rnorm(0,sigma) : 0.0;
    X = pow(K,(1-S))*pow(X,S)*exp(logeps);"
    ),
    delta.t=1
  ),
  paramnames=c("r","K","sigma","tau"),
  obsnames="Y",
  statenames="X"
) -> Gompertz


## ----params-------------------------------------------------------------------
Gompertz |>
  coef() |>
  parmat(4) -> p
p["X.0",] <- c(0.5,0.9,1.1,1.5)


## ----sim1---------------------------------------------------------------------
Gompertz |>
  simulate(params=p,format="data.frame") |>
  ggplot(aes(x=time,y=X,group=.id,color=.id))+
  geom_line()+
  guides(color="none")+
  theme_bw()+
  labs(title="Gompertz model",subtitle="stochastic simulations")


## ----pf1----------------------------------------------------------------------
Gompertz |>
  pfilter(Np=500) |>
  replicate(n=10) -> pf

logmeanexp(sapply(pf,logLik),se=TRUE,ess=TRUE)


## ----comparison,cache=TRUE----------------------------------------------------
system.time(simulate(gompertz,nsim=10000,format="arrays"))
system.time(simulate(Gompertz,nsim=10000,format="arrays"))
system.time(pfilter(gompertz,Np=10000))
system.time(pfilter(Gompertz,Np=10000))

