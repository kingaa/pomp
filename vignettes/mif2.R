## ----prelims,echo=FALSE,cache=FALSE--------------------------------------
library(ggplot2)
library(plyr)
library(reshape2)
library(magrittr)
theme_set(theme_bw())
library(pomp)
stopifnot(packageVersion("pomp")>="2.0.9.1")
options(
  keep.source=TRUE,
  stringsAsFactors=FALSE,
  encoding="UTF-8",
  scipen=5
)
set.seed(1332379783L)


## ----gompertz-init,results="hide"----------------------------------------
library(pomp)
gompertz() -> gomp
theta <- coef(gomp)
theta.true <- theta

## ----gompertz-sim,include=FALSE------------------------------------------
gomp %<>%
  window(start=1) %>%
  simulate(seed=340398091L)


## ----gompertz-mif2-1,results='hide'--------------------------------------
library(foreach)
library(doParallel)
registerDoParallel()

estpars <- c("r","sigma","tau")


## ----gompertz-mif2-2-eval,eval=TRUE,purl=TRUE,include=FALSE--------------
bake(file="gompertz-mif2.rds",{
  library(doRNG)
  registerDoRNG(525386942)
  
  foreach(i=1:10,.inorder=FALSE) %dopar% {
    theta.guess <- theta.true
    theta.guess[estpars] <- rlnorm(
      n=length(estpars),
      meanlog=log(theta.guess[estpars]),
      sdlog=1
    )
    gomp %>%
      mif2(
        Nmif=50,
        params=theta.guess,
        rw.sd=rw.sd(r=0.02,sigma=0.02,tau=0.05),
        cooling.fraction.50=0.95,
        Np=2000
      ) %>%
      continue(Nmif=50,cooling.fraction=0.8) %>%
      continue(Nmif=50,cooling.fraction=0.6) %>%
      continue(Nmif=50,cooling.fraction=0.2) -> m1
    ll <- replicate(n=10,logLik(pfilter(m1,Np=10000)))
    list(mif=m1,ll=logmeanexp(ll,se=TRUE))
  } -> mf
}) -> mf

## ----gompertz-mif-3------------------------------------------------------
lls <- sapply(mf,getElement,"ll")
best <- which.max(sapply(mf,getElement,"ll")[1,])
theta.mif <- coef(mf[[best]]$mif)

replicate(10,logLik(pfilter(gomp,params=theta.mif,Np=10000))) %>%
  logmeanexp(se=TRUE) -> pf.loglik.mif


## ----kf,include=FALSE----------------------------------------------------
kalman.filter <- function (Y, X0, r, K, sigma, tau) {
  ntimes <- length(Y)
  sigma.sq <- sigma^2
  tau.sq <- tau^2
  cond.loglik <- numeric(ntimes)
  filter.mean <- numeric(ntimes)
  pred.mean <- numeric(ntimes)
  pred.var <- numeric(ntimes)
  m <- log(X0)
  v <- 0
  S <- exp(-r)
  for (k in seq_len(ntimes)) {
    pred.mean[k] <- M <- (1-S)*log(K) + S*m
    pred.var[k] <- V <- S*v*S+sigma.sq
    q <- V+tau.sq
    r <- log(Y[k])-M
    cond.loglik[k] <- dnorm(x=log(Y[k]),mean=M,sd=sqrt(q),log=TRUE)-log(Y[k])
    q <- 1/V+1/tau.sq
    filter.mean[k] <- m <- (log(Y[k])/tau.sq+M/V)/q
    v <- 1/q
  }
  list(
    pred.mean=pred.mean,
    pred.var=pred.var,
    filter.mean=filter.mean,
    cond.loglik=cond.loglik,
    loglik=sum(cond.loglik)
  )
}

##' 'kalman' evaluates gompertz likelihood parameters X0 and K are fixed at 1.
##' Other parameters are taken from x (or, if not in x, from params).
kalman <- function (x, object, params) {
  Y <- obs(object)
  p <- params
  p[names(x)] <- x
  X0 <- 1
  r <- p["r"]
  K <- 1
  sigma <- p["sigma"]
  tau <- p["tau"]
  -kalman.filter(Y, X0, r, K, sigma, tau)$loglik
}

##' Exact log likelihood at the true parameters
loglik.truth <- -kalman(coef(gomp),gomp,coef(gomp))
loglik.mif <- -kalman(theta.mif,gomp,coef(gomp))

kalm.fit1 <- optim(
  par=theta.true[estpars],
  fn=kalman,
  object=gomp,
  params=coef(gomp),
  hessian=TRUE,
  control=list(trace=0)
)

theta.mle <- coef(gomp)
theta.mle[estpars] <- kalm.fit1$par
loglik.mle <- -kalm.fit1$value


## ----gomp-post-mif2,include=FALSE----------------------------------------
replicate(n=10,logLik(pfilter(gomp,Np=10000))) %>%
  logmeanexp(se=TRUE) -> pf.loglik.truth
replicate(n=10,logLik(pfilter(gomp,params=theta.mle,Np=10000))) %>%
  logmeanexp(se=TRUE) -> pf.loglik.mle

rbind(`Truth`=theta.true[estpars],
  `Exact MLE`=theta.mle[estpars],
  `IF2 MLE`=theta.mif[estpars]) %>% round(digits=4) %>%
  cbind(`$\\hat{\\ell}$`=round(c(pf.loglik.truth[1],pf.loglik.mle[1],pf.loglik.mif[1]),2),
    `s.e. $\\hat{\\ell}$`=round(c(pf.loglik.truth[2],pf.loglik.mle[2],pf.loglik.mif[2]),2),
    `$\\ell$`=round(c(loglik.truth,loglik.mle,loglik.mif),2)
  ) -> results.table


## ----mif2-plot,echo=FALSE,cache=FALSE,fig.height=6-----------------------
op <- par(mfrow=c(4,1),mar=c(3,3,0,0),mgp=c(2,1,0),bty='l')
loglik <- sapply(mf,function(x)conv.rec(x$mif,"loglik"))
r <- sapply(mf,function(x)conv.rec(x$mif,"r"))
sigma <- sapply(mf,function(x)conv.rec(x$mif,"sigma"))
tau <- sapply(mf,function(x)conv.rec(x$mif,"tau"))
matplot(loglik,type='l',lty=1,xlab="",ylab=expression(log~L),xaxt='n',ylim=max(loglik,na.rm=T)+c(-12,3))
matplot(r,type='l',lty=1,xlab="",ylab=expression(r),xaxt='n')
matplot(sigma,type='l',lty=1,xlab="",ylab=expression(sigma),xaxt='n')
matplot(tau,type='l',lty=1,xlab="MIF iteration",ylab=expression(tau))
par(op)


## ----first-mif-results-table,echo=FALSE,cache=FALSE----------------------
kable(results.table)

