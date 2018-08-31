options(digits=3)
png(filename="dacca-%02d.png",res=100)

library(pomp)

pompExample(dacca)

set.seed(1420306530L)

plot(dacca)
rinit(dacca)
coef(dacca)

stopifnot(all.equal(coef(dacca)[1:22],
  partrans(dacca,coef(dacca,transform=TRUE),dir="from")[1:22]))

plot(simulate(window(dacca,end=1893,seed=1420306530L)),yax.flip=TRUE)
pf <- freeze(pfilter(window(dacca,end=1893),Np=1000),seed=1420306530L)
plot(pf)

dev.off()
