library(pomp)

pompExample(ou2)
spy(ou2)

## fix some parameters
p <- c(
       alpha.1=0.9,alpha.2=0,alpha.3=0,alpha.4=0.99,
       sigma.1=1,sigma.2=0,sigma.3=2,
       tau=1,
       x1.0=50,x2.0=-50
       )

ou2.sim <- simulate(ou2,params=p,nsim=100,seed=32043858)

coef(ou2,c('x1.0','x2.0')) <- c(-50,50)

ou2.sim <- simulate(ou2)
x <- simulate(ou2,nsim=3,states=T)
y <- simulate(ou2,nsim=3,obs=T)
z <- simulate(ou2,nsim=3,obs=T,states=T)
window(ou2.sim,end=50) -> ignore
window(ou2.sim,end=150) -> ignore
time(ignore) <- seq(2,50,by=2)
try(coef(ignore,transform=TRUE) <- c(3,3))
coef(ignore) <- numeric(0)
try(coef(ignore) <- c(3,2))
coef(ignore) <- c(a=3,b=2)
coef(ignore,transform=TRUE) <- c(a=3,b=2)

try(simulate(ou2,times=c(3,2,1)))
try(simulate(ou2,times=c(1,2,2)))
try(simulate(ou2,times=c(-1,0,1)))
try(simulate(ou2,times=numeric()))

set.seed(577639485L)

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
  y <- obs(object)
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

pompExample(ou2)

# true coefficients
p.truth <- coef(ou2)
cat("coefficients at `truth'\n")
print(p.truth[c('alpha.1','alpha.4','x1.0','x2.0')],digits=4)
cat("Kalman filter log likelihood at truth\n")
print(loglik.truth <- -kalman(p.truth,ou2,p.truth),digits=4)

# make a wild guess
p.guess <- p.truth[c('alpha.1','alpha.4','x1.0','x2.0')]*exp(rnorm(n=4,mean=0,sd=0.5))
cat("coefficients at guess\n")
print(p.guess,digits=4)
cat("Kalman filter log likelihood at guess\n")
print(loglik.guess <- -kalman(p.guess,ou2,p.truth),digits=4)

# find MLE using Kalman filter starting at the guess
cat("running Kalman filter estimation\n")
kalm.fit1 <- optim(p.guess,kalman,object=ou2,params=p.truth,hessian=T,control=list(trace=2))
print(loglik.mle <- -kalm.fit1$value,digits=4)

cat("summary of results\n")
print(
      cbind(
            truth=c(p.truth[names(kalm.fit1$par)],loglik=loglik.truth),
            MLE=c(kalm.fit1$par,loglik=loglik.mle),
            SE=c(sqrt(diag(solve(kalm.fit1$hessian))),NA)
            ),
      digits=3
      )
