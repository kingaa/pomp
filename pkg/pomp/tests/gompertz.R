library(pomp)
options(digits=4)

pompExample(gompertz)

po <- gompertz
coef(po)
coef(po,transform=TRUE)
coef(po,c("r","X.0"))
coef(po,c("r","X.0"),transform=TRUE)
coef(po,c("r","K")) <- c(0.2,1)
coef(po)
coef(po,transform=TRUE)
guess <- coef(po)
guess["r"] <- 0

try(
    mf <- mif(
              po,
              start=guess,
              Nmif=5,Np=1000,
              transform=TRUE,
              ic.lag=1,var.factor=1,
              cooling.fraction=0.99^50,
              rw.sd=c(r=0.02,K=0.02)
              )
    )

set.seed(93848585L)
mf <- mif(
          po,
          Nmif=5,Np=1000,
          transform=TRUE,
          ic.lag=1,var.factor=1,
          cooling.fraction=0.99^50,
          rw.sd=c(r=0.02,K=0.02)
          )
coef(mf,transform=TRUE)
coef(mf)
conv.rec(mf)
conv.rec(mf,transform=TRUE)
conv.rec(mf,c("loglik","r"))
try(conv.rec(mf,c("loglik","r"),transform=FALSE))
try(conv.rec(mf,c("loglik","r"),transform=TRUE))
conv.rec(mf,c("loglik","r"),transform=TRUE)
conv.rec(mf,c("loglik"),transform=TRUE)
conv.rec(mf,c("K"),transform=TRUE)
