options(digits=3)
png(filename="gompertz-%02d.png",res=100)

library(pomp2)

pompExample(gompertz)

set.seed(1438408329L)

rinit(gompertz)
coef(gompertz)

stopifnot(all.equal(coef(gompertz),partrans(gompertz,coef(gompertz,transform=TRUE),dir="from")))
plot(simulate(gompertz,seed=1438408329L))
pf <- freeze(pfilter(gompertz,Np=1000),seed=1438408329L)
plot(pf)
tj <- trajectory(gompertz,params=c(K=1,r=0.1,sigma=0.1,tau=0.1,X.0=3))
plot(time(gompertz),tj[,,],type="l")

dev.off()
