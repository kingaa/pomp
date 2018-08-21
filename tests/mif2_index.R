options(digits=3)

library(pomp)
library(magrittr)

pompExample(gompertz)

set.seed(1481104436)

gompertz %>%
  mif2(Nmif=4,Np=1000,
    .indices=seq.int(1000),
    rw.sd=rw.sd(r=0.02,K=0.02,sigma=0.02),
    cooling.fraction=0.5) %>%
  slot("indices") %>%
  table()

set.seed(962724905)

gompertz %>%
  mif2(Nmif=4,Np=100,
    .indices=as.list(seq.int(100)),
    rw.sd=rw.sd(r=0.02,K=0.02,sigma=0.02),
    cooling.fraction=0.5) %>%
  slot("indices") %>%
  unlist() %>%
  table()

try(mif2(gompertz,Nmif=1,Np=100,
  .indices=1:5,rw.sd=rw.sd(r=0.02,K=0.02,sigma=0.02),
  cooling.fraction=0.5))
