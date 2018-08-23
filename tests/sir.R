library(pomp)

pompExample(sir,envir=NULL)[[1]] -> po

set.seed(48832734L)

stopifnot(all.equal(coef(po),
  partrans(po,coef(po,transform=TRUE),dir="from")))

po <- simulate(po)
stopifnot(
  round(mean(obs(po)))==434,
  round(mean(states(po,"S")))==133174)

pf <- pfilter(window(po,end=0.5),Np=1000)
stopifnot(
  round(mean(eff.sample.size(pf)))==558)

tj <- trajectory(po,maxsteps=10000,format="data.frame")
stopifnot(
  round(mean(tj$cases))==701,
  mean(tj$W)==0,
  round(mean(tj$I))==1397)
