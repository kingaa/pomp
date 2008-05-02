require(pomp)

data(ou2)

p.truth <- c(
             alpha.1=0.9,alpha.2=0,alpha.3=0,alpha.4=0.99,
             sigma.1=1,sigma.2=0,sigma.3=2,
             tau=1,
             x1.0=50,x2.0=-50
             )

ou2 <- mif(ou2,Nmif=0,start=p.truth,
           pars=c('alpha.1','alpha.4'),ivps=c('x1.0','x2.0'),
           rw.sd=c(
             x1.0=5,x2.0=5,
             alpha.1=0.1,alpha.2=0,alpha.3=0,alpha.4=0.1,
             sigma.1=0,sigma.2=0,sigma.3=0,
             tau=0
             ),
           alg.pars=list(Np=1000,var.factor=1,ic.lag=10,cooling.factor=0.95),
           max.fail=100
           )

fit1.pfilter <- pfilter(ou2)
cat("coefficients at `truth'\n")
print(coef(ou2,c('x1.0','x2.0','alpha.1','alpha.4')))
cat("particle filter log likelihood at truth\n")
print(fit1.pfilter$loglik)

guess <- ou2
coef(guess,c('x1.0','x2.0','alpha.1','alpha.4')) <- c(45,-60,0.8,0.9)

fit2.pfilter <- pfilter(guess,max.fail=1000,warn=F)
cat("coefficients at guess\n")
print(coef(guess,c('x1.0','x2.0','alpha.1','alpha.4')))
cat("particle filter log likelihood at guess\n")
print(fit2.pfilter$loglik)

cat("running MIF\n")
tic <- Sys.time()
mif.fit <- mif(guess,Nmif=80,max.fail=100)
toc <- Sys.time()
print(toc-tic)
cat("PF estimated log likelihood at MIF MLE\n")
print(pfilter(mif.fit,Np=10000)$loglik)

cat("coefficients at truth\n")
print(coef(ou2,c('x1.0','x2.0','alpha.1','alpha.4')))
cat("MIF MLE\n")
print(coef(mif.fit,c('x1.0','x2.0','alpha.1','alpha.4')))
