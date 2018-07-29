library(pomp)

pompExample(sir2,envir=NULL)[[1]] -> po
po <- window(po,end=1)

set.seed(48832734L)

stopifnot(all.equal(coef(po),partrans(po,coef(po,transform=TRUE),dir="from")))
po <- simulate(po)
stopifnot(sum(obs(po))==133)
stopifnot(round(mean(states(po,"S")),2)==56282.28)
pf <- pfilter(window(po,end=0.5),Np=1000)
stopifnot(all.equal(round(range(eff.sample.size(pf)),1),c(159.3,1000.0)))
tj <- trajectory(po,maxsteps=10000)
stopifnot(all.equal(round(apply(tj,1,mean),1),
  c(S=56327.1,I=38.7,R=943634.2,N=1e6,cases=17.2)))
