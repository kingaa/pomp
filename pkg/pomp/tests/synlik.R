library(pomp)

data(ou2)

set.seed(6457673L)

kalman.filter <- function (y, x0, a, b, sigma, tau) {
  n <- nrow(y)
  ntimes <- ncol(y)
  sigma.sq <- sigma%*%t(sigma)
  tau.sq <- tau%*%t(tau)
  inv.tau.sq <- solve(tau.sq)
  cond.dev <- numeric(ntimes)
  filter.mean <- matrix(0,n,ntimes)
  pred.mean <- matrix(0,n,ntimes)
  pred.var <- array(0,dim=c(n,n,ntimes))
  dev <- 0
  m <- x0
  v <- diag(0,n)
  for (k in seq(length=ntimes)) {
    pred.mean[,k] <- M <- a%*%m
    pred.var[,,k] <- V <- a%*%v%*%t(a)+sigma.sq
    q <- b%*%V%*%t(b)+tau.sq
    r <- y[,k]-b%*%M
    cond.dev[k] <- n*log(2*pi)+log(det(q))+t(r)%*%solve(q,r)
    dev <- dev+cond.dev[k]
    q <- t(b)%*%inv.tau.sq%*%b+solve(V)
    v <- solve(q)
    filter.mean[,k] <- m <- v%*%(t(b)%*%inv.tau.sq%*%y[,k]+solve(V,M))
  }
  list(
       pred.mean=pred.mean,
       pred.var=pred.var,
       filter.mean=filter.mean,
       cond.loglik=-0.5*cond.dev,
       loglik=-0.5*dev
       )
}

kalman <- function (x, object, params) {
  y <- data.array(object)
  p <- params
  p[names(x)] <- x
  x0 <- init.state(object,params=p)
  a <- matrix(p[c('alpha.1','alpha.2','alpha.3','alpha.4')],2,2)
  b <- diag(1,2)
  sigma <- matrix(p[c('sigma.1','sigma.2','sigma.2','sigma.3')],2,2)
  sigma[1,2] <- 0
  tau <- diag(p['tau'],2,2)
  -kalman.filter(y,x0,a,b,sigma,tau)$loglik
}

po <- window(ou2,end=5)

# exact likelihood
p.truth <- coef(po)
loglik.truth <- -kalman(p.truth,po,p.truth)

## likelihood from probes (works since ou2 is Gaussian)
loglik.probe <- replicate(n=500,logLik(probe(po,nsim=200,probes=function(x)x)))
## likelihood from particle filters
loglik.pfilter <- replicate(n=500,logLik(pfilter(po,Np=200)))

kruskal.test(list(loglik.probe,loglik.pfilter))
wilcox.test(loglik.probe,loglik.pfilter)
ks.test(loglik.pfilter,loglik.probe)
qqplot(loglik.pfilter,loglik.probe)
abline(a=0,b=1)
abline(v=loglik.truth)

##require(ggplot2)
##  x <- data.frame(probe=loglik.probe,pfilte=loglik.pfilter)
##  ggplot(data=melt(x))+geom_bar(aes(x=value,fill=variable),position="dodge")

