require(pomp.devel)
require(subplex)

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
  x0 <- p[c('x1.0','x2.0')]
  a <- matrix(p[c('alpha.1','alpha.2','alpha.3','alpha.4')],2,2)
  b <- diag(1,2)
  sigma <- matrix(p[c('sigma.1','sigma.2','sigma.2','sigma.3')],2,2)
  sigma[1,2] <- 0
  tau <- diag(p['tau'],2,2)
  -kalman.filter(y,x0,a,b,sigma,tau)$loglik
}

data(ou2)

# true coefficients
p.truth <- c(
             alpha.1=0.9,alpha.2=0,alpha.3=0,alpha.4=0.99,
             sigma.1=1,sigma.2=0,sigma.3=2,
             tau=1,
             x1.0=50,x2.0=-50
             )
cat("coefficients at `truth'\n")
print(p.truth[c('alpha.1','alpha.4','x1.0','x2.0')])
cat("Kalman filter log likelihood at truth\n")
print(-kalman(p.truth,ou2,p.truth))

# make a wild guess
p.guess <- c(alpha.1=0.8,alpha.4=0.9,x1.0=45,x2.0=-60)
cat("coefficients at guess\n")
print(p.guess)
cat("Kalman filter log likelihood at guess\n")
print(-kalman(p.guess,ou2,p.truth))

# find MLE using Kalman filter starting at the guess
cat("running Kalman filter estimation\n")
tic <- Sys.time()
kalm.fit1 <- optim(p.guess,kalman,object=ou2,params=p.truth,hessian=T)
toc <- Sys.time()
print(toc-tic)
tic <- Sys.time()
kalm.fit2 <- subplex(p.guess,function(x)kalman(x,ou2,p.truth),tol=1e-6)
toc <- Sys.time()
print(toc-tic)
cat("Kalman filter log likelihood at KF MLE\n")
print(-kalm.fit1$value)
print(-kalm.fit2$value)

cat("coefficients at truth\n")
print(p.truth[names(kalm.fit1$par)])
cat("Kalman filter MLE\n")
print(kalm.fit1$par)
print(kalm.fit2$par)
cat("Kalman filter estimated SEs\n")
print(sqrt(diag(solve(kalm.fit1$hessian))))
