options(digits=3)
png(filename="blowflies-%02d.png",res=100)

library(pomp)

capture.output(
  pompExample(blowflies,envir=NULL) -> flies,
  type="message") -> out
stopifnot(
  sum(grepl("unrecognized argument",out))==2,
  sum(grepl("y.init",out))==2
)

set.seed(599688L)

plot(flies[[1]])
rinit(flies[[1]])
coef(flies[[1]])
plot(simulate(flies[[1]],seed=599688L),var=c("y","R","S","N15"))
pf<- freeze(pfilter(flies[[1]],Np=1000),seed=599688L)
plot(pf)
stopifnot(
  partrans(
    flies[[1]],
    partrans(flies[[1]],dir="to",coef(flies[[1]])),
    dir="from"
  )==coef(flies[[1]]
  )
)

set.seed(1598688L)
plot(flies[[2]])
rinit(flies[[2]])
coef(flies[[2]])
plot(simulate(flies[[2]],seed=599688L),var=c("y","R","S","N8"))
pf <- freeze(pfilter(flies[[2]],Np=1000),seed=599688L)
plot(pf)
stopifnot(
  round(logLik(pf),2)==-1469.89,
  partrans(
    flies[[2]],
    partrans(flies[[2]],dir="to",coef(flies[[2]])),
    dir="from"
  )==coef(flies[[2]]
  )
)

dev.off()
