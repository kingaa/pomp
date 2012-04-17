require(pompExamples)

data(budmoth.sim)

names(budmoth.sim)
x <- lapply(budmoth.sim,as,"data.frame")

print(lapply(x,tail))

y <- simulate(budmoth.sim$food,seed=3434996L,as.data.frame=TRUE)
tail(y)

z <- trajectory(budmoth.sim$tri,as.data.frame=TRUE)
tail(z)

pf <- pfilter(budmoth.sim$food,seed=34348885L,Np=1000)
logLik(pf)
