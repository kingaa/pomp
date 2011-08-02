library(pomp)
options(digits=4)

data(gompertz)

po <- gompertz
coef(po)
coef(po,transform=TRUE)
coef(po,c("log.r","X.0"))
coef(po,c("r","X.0"),transform=TRUE)
coef(po,c("r","K"),transform=TRUE) <- c(0.2,1)
coef(po)
coef(po,transform=TRUE)

set.seed(93848585L)
mf <- mif(
          po,
          Nmif=5,Np=1000,
          ic.lag=1,var.factor=1,cooling.factor=0.99,
          rw.sd=c(log.r=0.02,log.K=0.02)
          )
coef(mf,transform=TRUE)
conv.rec(mf)
conv.rec(mf,transform=TRUE)
conv.rec(mf,c("loglik","log.r"))
try(conv.rec(mf,c("loglik","r"),transform=FALSE))
try(conv.rec(mf,c("loglik","log.r"),transform=TRUE))
conv.rec(mf,c("loglik","r"),transform=TRUE)
conv.rec(mf,c("loglik"),transform=TRUE)
conv.rec(mf,c("K"),transform=TRUE)

