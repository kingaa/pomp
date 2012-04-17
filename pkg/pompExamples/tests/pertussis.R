require(pompExamples)

data(pertussis.sim)

names(pertussis.sim)
x <- lapply(pertussis.sim,as.data.frame)

print(lapply(x,tail))

x <- simulate(pertussis.sim$full.big,seed=395885L,as.data.frame=TRUE)
tail(x)

y <- trajectory(pertussis.sim$SEIRS.small,as.data.frame=TRUE)
tail(y)

x <- simulate(pertussis.sim$full.small,seed=395885L,times=100,states=TRUE)[1:6,,1]
coef(pertussis.sim$full.small,paste0(names(x),".0")) <- x/sum(x)

pf <- pfilter(pertussis.sim$full.small,seed=3445886L,Np=1000)
logLik(pf)
