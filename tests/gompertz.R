options(digits=3)
png(filename="gompertz-%02d.png",res=100)

library(pomp)

gompertz() -> po

set.seed(1438408329L)

rinit(po)
coef(po)

stopifnot(all.equal(coef(po),partrans(po,coef(po,transform=TRUE),dir="from")))
plot(simulate(po,seed=1438408329L))
pf <- freeze(pfilter(po,Np=1000),seed=1438408329L)
plot(pf)
tj <- trajectory(po,params=c(K=1,r=0.1,sigma=0.1,tau=0.1,X_0=3),format="a")
plot(time(po),tj[,,],type="l")
e <- emeasure(po,x=states(po),params=coef(po),times=time(po))
plot(time(po),obs(po),xlab="",ylab="")
lines(time(po),t(apply(e,c(1,3),mean)))

dev.off()
