library(pomp)

## a stochastic version of the Verhulst-Pearl logistic model
## this evolves in continuous time, but we approximate the
## stochastic dynamics using an Euler approximation
## (plugin 'euler.sim') with fixed step-size

po <- pomp(
  data=data.frame(
    N=rep(0,1000),
    t=seq(0.1,by=0.1,length=1000)),
  times="t",
  t0=0,
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
  statenames=c("n")
)

params <- c(n_0=10000,K=10000,r=0.9,sigma=0.4,tau=0.1)
set.seed(73658676)
po <- simulate(po,params=params)

params <- cbind(
  c(n_0=100,K=10000,r=0.2,sigma=0.4,tau=0.1),
  c(n_0=1000,K=11000,r=0.1,sigma=0.4,tau=0.1)
)
traj <- trajectory(po,params=params,as.data.frame=TRUE)
traj <- reshape(traj,dir="wide",idvar="time",timevar="traj")
sim <- simulate(po,params=params,as.data.frame=TRUE,seed=34597368L)
sim <- reshape(sim,dir="wide",idvar="time",timevar="sim")
matplot(range(time(po)),range(c(traj[-1],sim[-1])),type='n',bty='l',lty=1,xlab="time",ylab="n",main="simulations vs trajectories")
matlines(traj$time,traj[-1],type='l',bty='l',lty=1,xlab="time",ylab="n")
matlines(sim$time,sim[c("n.1","n.2")],type='l',bty='l',lty=1,xlab="time",ylab="n")
