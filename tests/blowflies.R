png(filename="blowflies-%02d.png",res=100)
set.seed(599688L)

library(pomp)

list(blowflies1(),blowflies2()) -> flies

plot(flies[[1]])
rinit(flies[[1]])
coef(flies[[1]])
plot(simulate(flies[[1]],seed=599688L),var=c("y","R","S","N15"))
pf <- freeze(pfilter(flies[[1]],Np=1000),seed=599688L)
plot(pf)
logLik(pf)
stopifnot(
  all.equal(
    partrans(
      flies[[1]],
      partrans(flies[[1]],dir="to",coef(flies[[1]])),
      dir="from"
    ),
    coef(flies[[1]])
  ),
  abs(logLik(pf)+1467.82)<0.05
)

plot(flies[[2]])
rinit(flies[[2]])
coef(flies[[2]])
plot(simulate(flies[[2]],seed=599688L),var=c("y","R","S","N8"))
pf <- freeze(pfilter(flies[[2]],Np=1000),seed=599688L)
plot(pf)
logLik(pf)
stopifnot(
  all.equal(
    partrans(
      flies[[2]],
      partrans(flies[[2]],dir="to",coef(flies[[2]])),
      dir="from"
    ),
    coef(flies[[2]])
  ),
  abs(logLik(pf)+1473.66)<0.05
)

dev.off()
