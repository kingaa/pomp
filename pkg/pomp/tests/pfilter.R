library(pomp)

data(ou2)

pf <- pfilter(ou2,Np=1000,seed=343439L)
print(coef(ou2,c('x1.0','x2.0','alpha.1','alpha.4')),digits=4)
cat("particle filter log likelihood at truth\n")
print(pf$loglik,digits=4)

data(euler.sir)
pf <- pfilter(euler.sir,Np=200,seed=394343L)
print(pf$loglik,digits=4)
