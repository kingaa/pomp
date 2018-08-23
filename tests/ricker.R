library(pomp)

pompExample(ricker)

set.seed(48832734L)

stopifnot(
  all.equal(
    coef(ricker),
    partrans(ricker,coef(ricker,transform=TRUE),dir="from")))

po <- simulate(ricker)
stopifnot(sum(obs(po))==2040)

pf <- pfilter(ricker,Np=1000)
stopifnot(all.equal(round(mean(eff.sample.size(pf))),581))

tj <- trajectory(ricker,format="data.frame")
stopifnot(all.equal(round(mean(tj$N)),4))
