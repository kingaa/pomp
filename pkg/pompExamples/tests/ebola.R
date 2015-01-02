library(pompExamples)

set.seed(47575684L)

pompExample(ebola)
ebolaModel(country="Guinea") -> po
pf <- pfilter(simulate(po),Np=100)
tj <- trajectory(po)

ebolaModel(country="SierraLeone",na.rm=TRUE,type='cum') -> po
pf <- pfilter(simulate(po),Np=100)
tj <- trajectory(po)
dd <- simulate(po,as.data.frame=TRUE,obs=TRUE)
dd$week <- dd$time
po <- ebolaModel(data=subset(dd,select=c(week,cases,deaths)))
