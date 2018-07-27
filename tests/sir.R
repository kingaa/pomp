library(pomp)

pompExample(sir,envir=NULL)[[1]] -> po

set.seed(48832734L)

stopifnot(all.equal(coef(po),partrans(po,coef(po,transform=TRUE),dir="from")))
po <- simulate(po)
stopifnot(sum(obs(po))==90311)
stopifnot(round(mean(states(po,"S")),2)==133173.94)
pf <- pfilter(window(po,end=0.5),Np=1000)
stopifnot(all.equal(round(range(eff.sample.size(pf)),1),c(232.3,711.0)))
tj <- trajectory(po)
stopifnot(all.equal(round(sum(tj["cases",,]),2),145732.73))
