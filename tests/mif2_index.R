options(digits=3)

library(pomp)

gompertz() -> gompertz

set.seed(1481104436)

gompertz |>
  mif2(Nmif=4,Np=1000,
    .indices=seq.int(1000),
    rw.sd=rw_sd(r=0.02,K=0.02,sigma=0.02),
    cooling.fraction=0.5) |>
  slot("indices") -> idx
stopifnot(
  length(idx)==1000,
  class(idx)=="integer"
)

set.seed(962724905)

gompertz |>
  mif2(Nmif=4,Np=100,
    .indices=as.list(seq.int(100)),
    rw.sd=rw_sd(r=0.02,K=0.02,sigma=0.02),
    cooling.fraction=0.5) |>
  slot("indices") -> idx
stopifnot(
  length(idx)==100,
  class(idx)=="list"
)

try(mif2(gompertz,Nmif=1,Np=100,
  .indices=1:5,rw.sd=rw_sd(r=0.02,K=0.02,sigma=0.02),
  cooling.fraction=0.5))
