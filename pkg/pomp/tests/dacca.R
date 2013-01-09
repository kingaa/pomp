library(pomp)

set.seed(1420306530L)

pompExample(dacca)

x <- as.data.frame(dacca)
print(names(x))
print(dim(x))

x <- simulate(dacca,nsim=3,as.data.frame=TRUE)
print(names(x))
print(dim(x))

pf <- pfilter(dacca,Np=1000,seed=5886855L)
print(round(logLik(pf),digits=1))
