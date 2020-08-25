options(digits=3)
png(filename="ebola-%02d.png",res=100)

library(pomp)

ebolaModel() -> eb

set.seed(48832734L)

plot(eb)
coef(eb)
rinit(eb)

ebolaModel(country="SLE") -> eb

stopifnot(all.equal(coef(eb),
  partrans(eb,coef(eb,transform=TRUE),dir="from")))

plot(simulate(eb,seed=48832734L))
pf <- freeze(pfilter(window(eb,end=200),Np=1000),seed=48832734L)
plot(pf)
tj <- trajectory(eb,ode_control=list(maxsteps=10000))
matplot(time(eb),t(tj[c("I","N_EI","N_IR"),1,]),type="l",ylab="")

dev.off()
