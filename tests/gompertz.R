options(digits=3)

library(pomp)
library(magrittr)

pompExample(gompertz)

set.seed(1438408329L)

stopifnot(all.equal(coef(gompertz),partrans(gompertz,coef(gompertz,transform=TRUE),dir="from")))
po <- simulate(gompertz)
stopifnot(round(sum(obs(po)),3)==114.451)
pf <- pfilter(gompertz,Np=1000)
stopifnot(all.equal(round(mean(eff.sample.size(pf)),3),593.45))
tj <- trajectory(gompertz)
stopifnot(all.equal(round(sum(tj["X",,]),3),101))
