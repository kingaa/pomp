library(pomp)

pompExample(sir,envir=NULL)[[1]] -> po

set.seed(48832734L)

stopifnot(all.equal(coef(po),partrans(po,coef(po,transform=TRUE),dir="from")))
po <- simulate(po)
stopifnot(sum(obs(po))==90311)
stopifnot(round(mean(states(po,"S")),2)==133173.94)
pf <- pfilter(window(po,end=0.5),Np=1000)
stopifnot(all.equal(round(range(eff.sample.size(pf)),1),c(232.3,711.0)))
tj <- trajectory(po,maxsteps=10000)
stopifnot(all.equal(round(sum(tj["cases",,]),1),145732.7))
tj <- trajectory(po,maxsteps=10000)
stopifnot(all.equal(round(apply(tj,1,sum),1),
  c(S=27988531.2,I=290650.3,R=408520818.4,cases=145732.7,W=0)))
