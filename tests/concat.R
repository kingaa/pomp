options(digits=3)

library(pomp2)

try(pomp2:::concat())
try(pomp2:::concat("a","b"))

gompertz() -> gompertz
ou2() -> ou2
pomp2:::concat(a=ou2,c(b=gompertz,c=ou2))
c(a=ou2,c(b=gompertz,c=ou2)) %>% class()

replicate(2,pfilter(gompertz,Np=10)) %>% class()
do.call(c,replicate(2,pfilter(gompertz,Np=10))) -> pfs
pfs %>% class()
c(a=pfs[[1]],b=pfs)
