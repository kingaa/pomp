library(pomp)

data(ou2)

set.seed(64857673)

fit1.pfilter <- pfilter(ou2,Np=1000)
cat("coefficients at `truth'\n")
print(coef(ou2,c('x1.0','x2.0','alpha.1','alpha.4')),digits=4)
cat("particle filter log likelihood at truth\n")
print(fit1.pfilter$loglik,digits=4)

p.truth <- coef(ou2)
guess <- p.truth
guess[c('x1.0','x2.0','alpha.1','alpha.4')] <- c(45,-60,0.8,0.9)

fit2.pfilter <- pfilter(ou2,params=guess,Np=1000,max.fail=1000,warn=F)
cat("coefficients at guess\n")
print(guess[c('x1.0','x2.0','alpha.1','alpha.4')],digits=4)
cat("particle filter log likelihood at guess\n")
print(fit2.pfilter$loglik,digits=4)

cat("running MIF\n")
tic <- Sys.time()
mif.fit <- mif(ou2,Nmif=10,start=guess,
               pars=c('alpha.1','alpha.4'),ivps=c('x1.0','x2.0'),
               rw.sd=c(
                 x1.0=5,x2.0=5,
                 alpha.1=0.1,alpha.4=0.1
                 ),
               alg.pars=list(
                 Np=1000,
                 var.factor=1,
                 ic.lag=10,
                 cooling.factor=0.95
                 ),
               max.fail=100
               )
mif.fit <- continue(mif.fit,Nmif=70,max.fail=100)
toc <- Sys.time()
print(toc-tic)
cat("PF estimated log likelihood at MIF MLE\n")
print(pfilter(mif.fit)$loglik,digits=4)

cat("coefficients at truth\n")
print(coef(ou2,c('x1.0','x2.0','alpha.1','alpha.4')),digits=4)
cat("MIF MLE\n")
print(coef(mif.fit,c('x1.0','x2.0','alpha.1','alpha.4')),digits=4)

plot(mif.fit)
compare.mif(mif.fit)
compare.mif(list(mif.fit,mif.fit))

set.seed(33848585L)

fit <- mif(
           ou2,
           Nmif=0,
           pars=c("alpha.1","alpha.4"),
           ivps=c("x1.0","x2.0"),
           rw.sd=c(x1.0=5,x2.0=5,alpha.1=0.1,alpha.4=0.2,alpha.3=0),
           alg.pars=list(Np=100,cooling.factor=0.95,ic.lag=10,var.factor=1)
           )
fit <- mif(
           ou2,
           Nmif=1,
           pars=c("alpha.1","alpha.4"),
           ivps=c("x1.0","x2.0"),
           rw.sd=c(x1.0=5,x2.0=5,alpha.1=0.1,alpha.4=0.2,alpha.3=0),
           alg.pars=list(Np=1000,cooling.factor=0.95,ic.lag=10,var.factor=1)
           )
fit <- mif(
           ou2,
           Nmif=2,
           ivps=c("x1.0","x2.0"),
           rw.sd=c(x1.0=5,x2.0=5,alpha.1=0.1,alpha.4=0.2),
           alg.pars=list(Np=1000,cooling.factor=0.95,ic.lag=10,var.factor=1)
           )
fit <- continue(fit,Nmif=40)
ff <- pfilter(fit,pred.mean=T,filter.mean=T,pred.var=T,max.fail=100)
print(ff$loglik)
fit <- mif(fit,rw.sd=c(x1.0=5,x2.0=5,alpha.1=0.1,alpha.4=0.1))
fit <- continue(fit,Nmif=2,ivps=c("x1.0"),pars=c("alpha.1"))
s <- coef(fit)
s[2] <- 0.01
fit <- mif(fit,Nmif=10,start=s)
fit <- mif(ou2,Nmif=3,rw.sd=c(alpha.1=0.1,alpha.4=0.1),alg.pars=list(Np=1000,cooling.factor=0.98,var.factor=1,ic.lag=2))
fit <- continue(fit,Nmif=5,alg.pars=list(Np=2000,cooling.factor=0.98,var.factor=1,ic.lag=2))
fit <- continue(fit,ivps=c("x1.0"),rw.sd=c(alpha.1=0.1,alpha.4=0.1,x1.0=5,x2.0=5),Nmif=3)
