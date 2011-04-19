library(pomp)

data(ou2)

set.seed(6457673L)

po <- window(ou2,end=5)

pb <- replicate(n=100,logLik(probe(po,nsim=1000,probes=function(x)x)))
pf <- replicate(n=100,logLik(pfilter(po,Np=1000)))

kruskal.test(list(pb,pf))
ks.test(pf,pb)
qqplot(pf,pb)
abline(a=0,b=1)
