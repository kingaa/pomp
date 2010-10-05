library(pomp)

data(ricker)

set.seed(64857673L)

po <- ricker
coef(po) <- c(log.r=log(2),log.sigma=log(1),log.phi=log(20),N.0=7,e.0=0)

pf.tr <- pfilter(ricker,Np=1000,max.fail=50)
pf.po <- pfilter(po,Np=1000,max.fail=50)

mf <- mif(
          po,
          Nmif=100,
          Np=1000,
          cooling.factor=0.99,
          var.factor=2,
          ic.lag=3,
          max.fail=50,
          rw.sd=c(log.r=0.1,log.sigma=0.1,log.phi=0.1)
          )
plot(mf <- mif(mf,Nmif=500,max.fail=20))
pf.mf <- pfilter(mf,Np=1000)

res <- rbind(
             cbind(truth=coef(ricker),MLE=coef(mf),guess=coef(po)),
             loglik=c(pf.tr$loglik,pf.mf$loglik,pf.po$loglik)
             )

print(res,digits=3)

tj.1 <- trajectory(ricker)
plot(time(ricker),tj.1[1,,],type='l')
tj.2 <- trajectory(ricker,times=c(30:50),t0=0)
lines(30:50,tj.2[1,,],col='red',lwd=2)
max(abs(tj.1[,,time(ricker)>=30]-tj.2[,,]))

data(ricker)
po <- ricker
try(
    coef(po,"r")
    )
coef(po,c("r","phi")) <- c(0,0)
coef(po,c("log.r","log.phi")) <- c(a=0,b=0)
coef(po,c("log.r","log.phi")) <- 0
coef(po) <- c(log.phi=0,log.r=3.5,N.0=10,e.0=0,log.sigma=-Inf)
coef(po)
coef(po,"new") <- 3
plot(simulate(po))
coef(po)
