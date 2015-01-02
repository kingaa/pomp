library(pompExamples)

set.seed(47575684L)

pompExample(bbp)
pf <- pfilter(simulate(bbp),Np=100,max.fail=Inf)
tj <- trajectory(bbp)
