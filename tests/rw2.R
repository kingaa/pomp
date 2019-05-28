options(digits=3)
png(filename="rw2-%02d.png",res=100)

library(pomp)

rw2() -> rw2

set.seed(1438408329L)

rinit(rw2)
coef(rw2)

stopifnot(all.equal(coef(rw2),partrans(rw2,coef(rw2,transform=TRUE),dir="from")))
plot(simulate(rw2,seed=1438408329L))
pf <- freeze(pfilter(rw2,Np=1000),seed=1438408329L)
plot(pf)
tj <- trajectory(rw2)
try(matplot(time(rw2),t(tj[,1,]),type="l",ylab=""))

dev.off()
