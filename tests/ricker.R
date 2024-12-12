library(pomp)

set.seed(1800076828)

png(filename="ricker-%02d.png",res=100)

ricker() -> po
plot(po)

stopifnot(
  all.equal(
    coef(po),
    partrans(po,coef(po,transform=TRUE),dir="from")
  ),
  all.equal(
    coef(po,transform=TRUE),
    partrans(po,coef(po),dir="to")
  )
)

pfilter(
  po,
  Np=1000,
  filter.mean=TRUE,
  pred.mean=TRUE,
  pred.var=TRUE,
  filter.traj=TRUE,
  save.states="filter"
) -> pf

stopifnot(
  abs(logLik(pf)+138.5)<0.05
)

plot(pf,yax.flip=TRUE)

forecast(pf,format="d") -> fc
simulate(pf) -> sm

emeasure(sm) -> ef
vmeasure(sm) -> vf
plot(ef,vf)

plot(time(sm),obs(sm),xlab="time",ylab="Y")
lines(time(sm),ef)

trajectory(po) -> tj
plot(tj)

dev.off()
