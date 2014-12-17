library(pompExamples)

## pdf.options(useDingbats=FALSE)
## pdf(file="examples.pdf")

set.seed(47575684L)

po <- pompExample(parus,proc="Ricker",meas="lognormal",envir=NULL)
pf <- pfilter(simulate(po$parus),Np=100,max.fail=Inf)
tj <- trajectory(po$parus)

po <- pompExample(parus,proc="Ricker",meas="negbin",envir=NULL)
pf <- pfilter(simulate(po$parus),Np=100,max.fail=Inf)

po <- pompExample(parus,proc="Ricker",meas="Poisson",envir=NULL)
pf <- pfilter(simulate(po$parus),Np=100,max.fail=Inf)

po <- pompExample(parus,proc="Gompertz",meas="Poisson",envir=NULL)
pf <- pfilter(simulate(po[[1]]),Np=100,max.fail=Inf)
tj <- trajectory(po[[1]])

po <- pompExample(parus,proc="Gompertz",meas="lognormal",envir=NULL)
pf <- pfilter(simulate(po$parus),Np=100,max.fail=Inf)

po <- pompExample(bbp,envir=NULL)
pf <- pfilter(simulate(pf$bbp),Np=100,max.fail=Inf)

## dev.off()
