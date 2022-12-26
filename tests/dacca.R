library(pomp)

set.seed(1420306530L)

png(filename="dacca-%02d.png",res=100)

dacca() -> po
plot(po,yax.flip=TRUE)
coef(po)

theta1 <- coef(po)
theta2 <- partrans(po,coef(po,transform=TRUE),dir="from")

stopifnot(
  all.equal(
    theta1[1:22],
    theta2[1:22]
  ),
  all.equal(
    theta1[23:28]/sum(theta1[23:28]),
    theta2[23:28]/sum(theta2[23:28])
  ),
  all.equal(
    coef(po,transform=TRUE),
    partrans(po,coef(po),dir="to")
  )
)

pfilter(
  window(po,end=1893),
  Np=1000,
  filter.mean=TRUE,
  pred.mean=TRUE,
  pred.var=TRUE,
  filter.traj=TRUE,
  save.states=TRUE
) -> pf

stopifnot(
  abs(logLik(pf)+150.3)<0.5
)

plot(pf,yax.flip=TRUE)

pf %>% window(end=1893) %>% simulate() -> sm
plot(sm,yax.flip=TRUE)

try(dacca(logbeta=c(1,2,3),logomega=c(10,20)))

dev.off()
