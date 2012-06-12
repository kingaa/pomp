library(pomp)

data(ou2)

set.seed(9994847L)

pf <- pfilter(ou2,Np=1000,seed=343439L)
print(coef(ou2,c('x1.0','x2.0','alpha.1','alpha.4')),digits=4)
cat("particle filter log likelihood at truth\n")
print(pf$loglik,digits=4)

pf <- replicate(n=10,pfilter(ou2,Np=function(k)if(k<10) 10000 else 500))
pf.ll <- sapply(pf,logLik)
ll.est <- log(mean(exp(pf.ll-mean(pf.ll))))+mean(pf.ll)
ll.se <- sd(exp(pf.ll-mean(pf.ll)))/exp(ll.est-mean(pf.ll))/sqrt(length(pf))
print(round(c(loglik=ll.est,loglik.se=ll.se),digits=2))

data(euler.sir)
pf <- pfilter(euler.sir,Np=100,seed=394343L)
print(coef(pf))
print(pf$loglik,digits=4)

p <- coef(euler.sir)
euler.sir@params <- numeric(0)
p["iota"] <- 1
pf <- pfilter(euler.sir,params=p,Np=100,seed=394343L)
print(coef(pf))
print(logLik(pf),digits=4)
