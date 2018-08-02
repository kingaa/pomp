options(digits=3)
library(pomp)

pompExample(dacca)

set.seed(1420306530L)

stopifnot(all.equal(coef(dacca)[1:22],
  partrans(dacca,coef(dacca,transform=TRUE),dir="from")[1:22]))

x <- simulate(window(dacca,end=1893),nsim=3,as.data.frame=TRUE)
stopifnot(all.equal(round(mean(x$cholera.deaths)),1003))

pf <- pfilter(window(dacca,end=1893),Np=1000)
stopifnot(abs(-150.9-logLik(pf))<2)
