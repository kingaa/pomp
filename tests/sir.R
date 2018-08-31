options(digits=3)
png(filename="sir-%02d.png",res=100)

library(pomp)

pompExample(sir)

set.seed(48832734L)

plot(sir)
coef(sir)
rinit(sir)

stopifnot(all.equal(coef(sir),
  partrans(sir,coef(sir,transform=TRUE),dir="from")))

plot(simulate(sir,seed=48832734L))
pf <- freeze(pfilter(window(sir,end=0.5),Np=1000),seed=48832734L)
plot(pf)
tj <- trajectory(sir,maxsteps=10000)
matplot(time(sir),t(tj[c("I","cases"),1,]),type="l",ylab="")

dev.off()
