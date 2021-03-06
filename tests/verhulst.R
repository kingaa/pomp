options(digits=3)
png(filename="verhulst-%02d.png",res=100)

library(pomp)

verhulst(n_0=100) -> po

set.seed(1438408329L)

rinit(po)
coef(po)
plot(po)

stopifnot(all.equal(coef(po),partrans(po,coef(po,transform=TRUE),dir="from")))
plot(simulate(po,seed=1438408329L))
pf <- freeze(pfilter(po,Np=1000),seed=1438408329L)
plot(pf)
tj <- trajectory(po)
plot(time(po),tj["n",1,],type="l",ylab="n")

dev.off()
