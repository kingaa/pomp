library(pomp)

set.seed(1438408329L)

png(filename="verhulst-%02d.png",res=100)

verhulst(n_0=100) -> po
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
  save.states='filter'
) -> pf

print(logLik(pf))
stopifnot(
  abs(logLik(pf)+8687.3)<0.05
)

plot(pf,yax.flip=TRUE)

forecast(pf,format="d") -> fc
simulate(pf) -> sm

emeasure(sm) -> ef
vmeasure(sm) -> vf
plot(ef,vf)

plot(time(sm),obs(sm),xlab="time",ylab="Y")
lines(time(sm),ef)

enkf(po,Np=1000) -> kf
plot(kf,yax.flip=TRUE)

trajectory(po) -> tj
plot(tj)

dev.off()
