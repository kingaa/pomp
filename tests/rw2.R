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
tj <- trajectory(rw2,format="array")
stopifnot(
  is.na(tj)
)

v <- vmeasure(rw2,params=coef(rw2),x=states(rw2),times=time(rw2))
stopifnot(
  v[1,1,,]==1,
  v[2,2,,]==1,
  v[1,2,,]==0,
  v[2,1,,]==0
)

dev.off()
