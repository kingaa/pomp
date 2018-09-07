library(pomp2)

## A stochastic version of the Verhulst-Pearl logistic model.
## This evolves in continuous time, but we approximate the
## stochastic dynamics using an Euler approximation
## (plugin 'euler') with fixed step-size ('delta.t').

simulate(
  times=seq(0.1,by=0.1,length=1000),
  t0=0,
  params=c(n_0=10000,K=10000,r=0.9,sigma=0.4,tau=0.1),
  rprocess=euler(
    step.fun=Csnippet("
        n = rnorm(n+r*n*(1-n/K)*dt,sigma*n*sqrt(dt));
      "
    ),
    delta.t=0.01
  ),
  dprocess=Csnippet("
        double dt = t_2-t_1;
        loglik = dnorm(n_2,n_1+r*n_1*(1-n_1/K)*dt,sigma*n_1*sqrt(dt),1);
      "
  ),
  rmeasure=Csnippet("N = rlnorm(log(n),log(1+tau));"),
  dmeasure=Csnippet("lik = dlnorm(N,log(n),log(1+tau),give_log);"),
  skeleton=vectorfield(Csnippet("Dn = r*n*(1-n/K);")),
  rinit=Csnippet("n = n_0;"),
  paramnames=c("r","K","tau","sigma","n_0"),
  statenames=c("n"),
  obsnames="N",
  seed=73658676L
) -> logistic

c("logistic")
