library(pomp)
library(magrittr)

## a simple two-dimensional random walk
## this makes use of the 'onestep.sim' plugin
## which we can use since we can simulate the
## increment of a random walk over any time

simulate(
  times=seq(1,100), t0=0,
  params=c(s1=3,s2=1,x1.0=0,x2.0=1,tau=10),
  rprocess = onestep.sim(
    step.fun = Csnippet("
      x1 = rnorm(x1,s1*sqrt(dt));
      x2 = rnorm(x2,s2*sqrt(dt));
      "
    )
  ),
  dprocess = Csnippet("
      double sdt = sqrt(t_2 - t_1);
      loglik = dnorm(x1_2,x1_1,s1*sdt,1)+
          dnorm(x2_2,x2_1,s2*sdt,1);
      "
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
  obsnames=c("y1","y2")
) -> rw2

plot(rw2)
pf <- pfilter(rw2,Np=1000)
logLik(pf)
plot(pf)
