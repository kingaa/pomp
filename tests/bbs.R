library(pomp)

pompExample(bbs)

set.seed(48832734L)

stopifnot(all.equal(coef(bbs),partrans(bbs,coef(bbs,transform=TRUE),dir="from")))
po <- simulate(bbs)
stopifnot(sum(obs(po))==583)
pf <- pfilter(bbs,Np=1000)
stopifnot(all.equal(round(mean(eff.sample.size(pf)),3),904.143))
tj <- trajectory(bbs)
stopifnot(all.equal(round(sum(tj["H",,]),3),1209.546))
