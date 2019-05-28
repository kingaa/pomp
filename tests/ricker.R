options(digits=3)
png(filename="ricker-%02d.png",res=100)

library(pomp)

ricker() -> ricker

set.seed(1438408329L)

rinit(ricker)
coef(ricker)
plot(ricker)

stopifnot(all.equal(coef(ricker),partrans(ricker,coef(ricker,transform=TRUE),dir="from")))
plot(simulate(ricker,seed=1438408329L))
pf <- freeze(pfilter(ricker,Np=1000),seed=1438408329L)
plot(pf)
tj <- trajectory(ricker)
plot(time(ricker),tj["N",1,],type="l",ylab="N")

dev.off()
