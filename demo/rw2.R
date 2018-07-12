library(pomp)

## a simple two-dimensional random walk
## this makes use of the 'onestep.sim' plugin
## which we can use since we can simulate the
## increment of a random walk over any time

rw2 <- pomp(
  rprocess = onestep.sim(
    step.fun = Csnippet("
      x1 = rnorm(x1,s1*dt);
      x2 = rnorm(x2,s2*dt);
      "
    )
  ),
  dprocess = onestep.dens(
    dens.fun = Csnippet("
      double dt = t_2 - t_1;
      loglik = dnorm(x1_2,x1_1,s1*dt,1)+
          dnorm(x2_2,x2_1,s2*dt,1);
      "
    )
  ),
  rmeasure=Csnippet("
      y1 = rnorm(x1,tau);
      y2 = rnorm(x2,tau);
    "
  ),
  dmeasure=Csnippet("
      lik = dnorm(y1,x1,tau,1)+dnorm(y2,x2,tau,1);
      lik = (give_log) ? lik : exp(lik);
    "
  ),
  statenames=c("x1","x2"),
  paramnames=c("s1","s2","tau"),
  data=data.frame(
    t=1:100,
    y1=rep(0,100),
    y2=rep(0,100)
  ),
  times='t',
  t0=0
)

rw2 <- simulate(rw2,params=c(s1=3,s2=1,x1.0=0,x2.0=1,tau=10))
plot(rw2)
pf <- pfilter(rw2,Np=100)
logLik(pf)
plot(pf)
