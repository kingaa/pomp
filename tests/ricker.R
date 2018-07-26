library(pomp)

pompExample(ricker)

set.seed(48832734L)

stopifnot(all.equal(coef(ricker),partrans(ricker,coef(ricker,transform=TRUE),dir="from")))
po <- simulate(ricker)
stopifnot(sum(obs(po))==2040)
pf <- pfilter(ricker,Np=1000)
stopifnot(all.equal(round(mean(eff.sample.size(pf)),3),580.979))
tj <- trajectory(ricker)
stopifnot(all.equal(round(sum(tj["N",,]),3),193.225))
