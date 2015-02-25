library(pompExamples)

set.seed(47575684L)

pompExample(parus)

po <- parus(proc="Ricker",meas="lognormal")
pf <- pfilter(simulate(po),Np=100,max.fail=Inf)
tj <- trajectory(po)

po <- parus(proc="Ricker",meas="negbin")
pf <- pfilter(simulate(po),Np=100,max.fail=Inf)

po <- parus(proc="Ricker",meas="Poisson")
pf <- pfilter(simulate(po),Np=100,max.fail=Inf)

po <- parus(proc="Gompertz",meas="Poisson")
pf <- pfilter(simulate(po),Np=100,max.fail=Inf)
tj <- trajectory(po)

po <- parus(proc="Gompertz",meas="lognormal")
pf <- pfilter(simulate(po),Np=100,max.fail=Inf)
