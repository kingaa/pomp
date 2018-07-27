options(digits=3)
library(pomp)

pompExample(dacca)

set.seed(1420306530L)

x <- simulate(dacca,nsim=3,as.data.frame=TRUE)
stopifnot(all.equal(round(sum(x$cholera.deaths),2),1011741.66))
pf <- pfilter(dacca,Np=1000)
stopifnot(all.equal(round(logLik(pf),2),-3753.56))
