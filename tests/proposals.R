library(pomp)
set.seed(1178744046L)

pompExample(ou2)

f <- mvn.diag.rw(c(a=10,10))
try(pmcmc(ou2,Nmcmc=2,Np=100,proposal=f))
try(abc(ou2,Nmcmc=2,Np=100,proposal=f,probes=list(probe.mean("y1")),
        scale=1,epsilon=1))
