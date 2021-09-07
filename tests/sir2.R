options(digits=3)
png(filename="sir2-%02d.png",res=100)

library(pomp)

sir2() -> sir2

set.seed(48832734L)

plot(sir2)
coef(sir2)
rinit(sir2)

stopifnot(all.equal(coef(sir2),
  partrans(sir2,coef(sir2,transform=TRUE),dir="from")))

plot(simulate(sir2,seed=48832734L))
pf <- freeze(pfilter(window(sir2,end=0.5),Np=1000),seed=48832734L)
plot(pf)
tj <- trajectory(sir2,ode_control=list(maxsteps=10000),format="a")
matplot(time(sir2),t(tj[c("I","cases"),1,]),type="l",ylab="")

dev.off()
