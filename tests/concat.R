options(digits=3)

library(pomp)

try(pomp:::concat())
try(pomp:::concat("a","b"))

gompertz() -> gompertz
ou2() -> ou2
pomp:::concat(a=ou2,c(b=gompertz,c=ou2))
c(a=ou2,c(b=gompertz,c=ou2)) |> class()
c(a=ou2,b=ou2(alpha_1=-11)) -> pomps
pomps |> coef()
pomps |> obs(vars="y1") |> melt() |> head()
pomps |> states(vars="x2") |> melt() |> head()

replicate(2,pfilter(gompertz,Np=10)) |> class()
do.call(c,replicate(2,pfilter(gompertz,Np=10))) -> pfs
pfs |> class()
c(a=pfs[[1]],b=pfs) -> pfs
pfs
time(pfs) -> tt
states(pfs) -> ss
obs(pfs) -> oo
coef(pfs) -> cc
stopifnot(
  is.list(tt),
  names(tt)==names(pfs),
  sapply(tt,length)==100,
  is.list(ss),
  names(ss)==names(pfs),
  sapply(ss,dim)==c(1,100),
  is.list(oo),
  names(oo)==names(pfs),
  sapply(oo,dim)==c(1,100),
  is.matrix(cc),
  dim(cc)==c(5,3),
  rownames(cc)==names(coef(pfs[[1]])),
  colnames(cc)==names(pfs)
)

