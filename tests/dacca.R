options(digits=3)
png(filename="dacca-%02d.png",res=100)

library(pomp)

dacca() -> po

set.seed(1420306530L)

plot(po)
rinit(po)
coef(po)

stopifnot(all.equal(coef(po)[1:22],
  partrans(po,coef(po,transform=TRUE),dir="from")[1:22]))

plot(simulate(window(po,end=1893,seed=1420306530L)),yax.flip=TRUE)
pf <- freeze(pfilter(window(po,end=1893),Np=1000),seed=1420306530L)
plot(pf)

try(dacca(logbeta=c(1,2,3),logomega=c(10,20)))

dev.off()
