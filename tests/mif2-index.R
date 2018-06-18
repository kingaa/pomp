library(pomp)
options(digits=4)

pompExample(gompertz)

set.seed(1481104436)

mf <- mif2(gompertz,Nmif=4,Np=1000,
           transform=TRUE,
           .indices=seq.int(1000),
           rw.sd=rw.sd(r=0.02,K=0.02,sigma=0.02),
           cooling.fraction=0.5)
table(mf@indices)

set.seed(962724905)

mf <- mif2(gompertz,Nmif=4,Np=100,
           transform=TRUE,
           .indices=as.list(seq.int(100)),
           rw.sd=rw.sd(r=0.02,K=0.02,sigma=0.02),
           cooling.fraction=0.5)
table(unlist(mf@indices))

try(mif2(gompertz,Nmif=1,Np=100,transform=TRUE,
  .indices=1:5,rw.sd=rw.sd(r=0.02,K=0.02,sigma=0.02),
  cooling.fraction=0.5))
