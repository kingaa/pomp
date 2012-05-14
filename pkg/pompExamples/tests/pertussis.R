require(pompExamples)

data(pertussis.sim)

names(pertussis.sim)
x <- lapply(pertussis.sim,as.data.frame)

print(lapply(x,tail))

x <- simulate(pertussis.sim$full.big,seed=395885L,as.data.frame=TRUE)
tail(x)

y <- trajectory(pertussis.sim$SEIRS.small,as.data.frame=TRUE)
tail(y)

system.time(pf <- pfilter(pertussis.sim$full.small,seed=3445886L,Np=1000))
logLik(pf)

pttest <- function (po, digits = 15) {
  identical(
            signif(coef(po),digits=digits),
            signif(partrans(po,partrans(po,coef(po),dir='inv'),dir='for'),digits=digits)
            )
}

stopifnot(all(sapply(pertussis.sim,pttest)))

pttest <- function (po, digits = 15) {
  identical(
            signif(coef(po,trans=T),digits=digits),
            signif(partrans(po,partrans(po,coef(po,trans=T),dir='f'),dir='inv'),digits=digits)
            )
}

stopifnot(all(sapply(pertussis.sim,pttest)))

