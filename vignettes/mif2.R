## ----prelims,echo=FALSE,cache=FALSE--------------------------------------
require(ggplot2)
require(plyr)
require(reshape2)
require(magrittr)
theme_set(theme_bw())
require(pomp)
stopifnot(packageVersion("pomp")>="0.66-2")
options(
  keep.source=TRUE,
  stringsAsFactors=FALSE,
  encoding="UTF-8",
  scipen=5,
  cores=5
  )

## ----gompertz-init,cache=FALSE-------------------------------------------
require(pomp)
pompExample(gompertz)
theta <- coef(gompertz)
theta.true <- theta

## ----gompertz-multi-mif2-eval,results='hide'-----------------------------
require(foreach)
require(doMC)
registerDoMC()

save.seed <- .Random.seed
set.seed(334388458L,kind="L'Ecuyer")

estpars <- c("r","sigma","tau")
mf <- foreach(i=1:10,
              .inorder=FALSE,
              .options.multicore=list(set.seed=TRUE)
              ) %dopar%
  {
    theta.guess <- theta.true
    theta.guess[estpars] <- rlnorm(
      n=length(estpars),
      meanlog=log(theta.guess[estpars]),
      sdlog=1
      )
    m1 <- mif2(
      gompertz,
      Nmif=50,
      start=theta.guess,
      transform=TRUE,
      rw.sd=rw.sd(r=0.02,sigma=0.02,tau=0.05),
      cooling.fraction.50=0.95,
      Np=2000
      )
    m1 <- continue(m1,Nmif=50,cooling.fraction=0.8)
    m1 <- continue(m1,Nmif=50,cooling.fraction=0.6)
    m1 <- continue(m1,Nmif=50,cooling.fraction=0.2)
    ll <- replicate(n=10,logLik(pfilter(m1,Np=10000)))
    list(mif=m1,ll=ll)
    }

## ----gompertz-post-mif2--------------------------------------------------
loglik.true <- replicate(n=10,logLik(pfilter(gompertz,Np=10000)))
loglik.true <- logmeanexp(loglik.true,se=TRUE)
theta.mif <- t(sapply(mf,function(x)coef(x$mif)))
loglik.mif <- t(sapply(mf,function(x)logmeanexp(x$ll,se=TRUE)))
best <- which.max(loglik.mif[,1])
theta.mif <- theta.mif[best,]
loglik.mif <- loglik.mif[best,]
rbind(
  mle=c(signif(theta.mif[estpars],3),loglik=round(loglik.mif,2)),
  truth=c(signif(theta.true[estpars],3),loglik=round(loglik.true,2))
  ) -> results.table

## ----mif2-plot,echo=FALSE,cache=FALSE,fig.height=6-----------------------
op <- par(mfrow=c(4,1),mar=c(3,3,0,0),mgp=c(2,1,0),bty='l')
loglik <- sapply(mf,function(x)conv.rec(x$mif,"loglik"))
log.r <- sapply(mf,function(x)conv.rec(x$mif,"r"))
log.sigma <- sapply(mf,function(x)conv.rec(x$mif,"sigma"))
log.tau <- sapply(mf,function(x)conv.rec(x$mif,"tau"))
matplot(loglik,type='l',lty=1,xlab="",ylab=expression(log~L),xaxt='n',ylim=max(loglik,na.rm=T)+c(-12,3))
matplot(log.r,type='l',lty=1,xlab="",ylab=expression(log~r),xaxt='n')
matplot(log.sigma,type='l',lty=1,xlab="",ylab=expression(log~sigma),xaxt='n')
matplot(log.tau,type='l',lty=1,xlab="MIF iteration",ylab=expression(log~tau))
par(op)

## ----first-mif-results-table,echo=FALSE,cache=FALSE----------------------
print(results.table)

