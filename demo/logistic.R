library(pomp)
library(magrittr)
library(dplyr)
library(ggplot2)

## a stochastic version of the Verhulst-Pearl logistic model
## this evolves in continuous time, but we approximate the
## stochastic dynamics using an Euler approximation
## (plugin 'euler.sim') with fixed step-size

simulate(
  times=seq(0.1,by=0.1,length=1000),
  t0=0,
  params=c(n_0=10000,K=10000,r=0.9,sigma=0.4,tau=0.1),
  rprocess=euler.sim(
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
) -> po

params <- cbind(
  c(n_0=100,K=10000,r=0.2,sigma=0.4,tau=0.1),
  c(n_0=100,K=20000,r=0.1,sigma=0.4,tau=0.1),
  c(n_0=100,K=30000,r=0.3,sigma=0.2,tau=0.1)
)

po %>% trajectory(params=params,as.data.frame=TRUE) %>%
  rename(rep=.id) %>%
  mutate(type="deterministic") -> traj
po %>% simulate(params=params,format="data.frame",seed=34597368L) %>%
  rename(rep=.id) %>%
  mutate(type="stochastic") -> sim
full_join(traj,sim) %>%
  ggplot(aes(x=time,y=n,group=interaction(rep,type),color=type))+
  geom_line()+
  theme_bw()+
  labs(title="Logistic model",subtitle="simulations vs. trajectories",color="")
