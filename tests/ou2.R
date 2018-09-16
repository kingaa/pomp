options(digits=3)
png(filename="ou2-%02d.png",res=100)

library(pomp2)

ou2() -> ou2

set.seed(1438408329L)

plot(ou2)
rinit(ou2)
coef(ou2)

stopifnot(all.equal(coef(ou2),partrans(ou2,coef(ou2,transform=TRUE),dir="from")))
plot(s <- simulate(ou2,seed=1438408329L))
pf <- freeze(pfilter(ou2,Np=1000),seed=1438408329L)
plot(pf)
tj <- trajectory(ou2)
matplot(time(ou2),t(tj[,1,]),type="l",ylab="")

d <- dprocess(s,x=states(s),params=coef(s),times=time(s),log=TRUE)
plot(d[,],ylab="log prob")

try(dprocess(s,x=states(s)[,c(1:9,15)],params=coef(s),times=time(s)[c(1:9,15)]))

dev.off()
