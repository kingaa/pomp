library(pompExamples)

all <- c("SEIR.small","SEIR.big","SEIRS.small","SEIRS.big","SEIRR.small","SEIRR.big","full.small","full.big")

sapply(all,function(n)eval(bquote(pertussis.sim(.(n))))) -> pt

names(pt)
x <- lapply(pt,as.data.frame)

print(lapply(x,tail))

x <- simulate(pertussis.sim(full.big),seed=395885L,as.data.frame=TRUE)
tail(x)

y <- trajectory(pertussis.sim(SEIRS.small),as.data.frame=TRUE)
tail(y)

pf <- pfilter(pertussis.sim(full.small),seed=3445886L,Np=1000)
logLik(pf)

pttest <- function (po, digits = 15) {
  identical(
            signif(coef(po),digits=digits),
            signif(partrans(po,partrans(po,coef(po),dir='inv'),dir='for'),digits=digits)
            )
}

stopifnot(all(sapply(pt,pttest)))

pttest <- function (po, digits = 15) {
  identical(
            signif(coef(po,trans=T),digits=digits),
            signif(partrans(po,partrans(po,coef(po,trans=T),dir='f'),dir='inv'),digits=digits)
            )
}

stopifnot(all(sapply(pt,pttest)))
