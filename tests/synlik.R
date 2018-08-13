# options(digits=3)
# png(filename="synlik-%02d.png",res=100)

library(pomp)

pompExample(ou2)

set.seed(6457673L)

po <- window(ou2,end=5)

ou2.kalman <- function (x, object, params) {
  y <- obs(object)
  p <- params
  p[names(x)] <- x
  x0 <- rinit(object,params=p)
  A <- matrix(p[c('alpha.1','alpha.2','alpha.3','alpha.4')],2,2)
  C <- diag(1,2)
  Q <- matrix(p[c('sigma.1','sigma.2',NA,'sigma.3')],2,2)
  Q[1,2] <- 0
  Q <- tcrossprod(Q)
  R <- diag(p['tau']^2,2,2)
  pomp:::kalmanFilter(t=time(object),y=y,X0=x0,A=A,C=C,Q=Q,R=R)$loglik
}

# exact likelihood
p.truth <- coef(po)
loglik.truth <- ou2.kalman(p.truth,po,p.truth)

## likelihood from probes (works since ou2 is Gaussian)
loglik.probe <- replicate(n=500,logLik(probe(po,nsim=200,probes=function(x)x)))
## likelihood from particle filters
loglik.pfilter <- replicate(n=500,logLik(pfilter(po,Np=200)))

kruskal.test(list(loglik.probe,loglik.pfilter))
wilcox.test(loglik.probe,loglik.pfilter)
ks.test(loglik.pfilter,loglik.probe)
plot(density(loglik.probe))
abline(v=loglik.truth)
plot(density(loglik.pfilter))
abline(v=loglik.truth)

# dev.off()
