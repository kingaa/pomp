library(pomp)
suppressPackageStartupMessages({
  library(dplyr)
})

png(filename="ebola-%02d.png",res=100)

set.seed(48832734L)

ebolaModel() -> po
plot(po)
coef(po)

ebolaWA2014 |>
  filter(
    country=="SLE",
    date<="2014-10-31"
  ) |>
  mutate(day=as.numeric(date-as.Date("2014-04-30"))) |>
  select(-date,-country) |>
  ebolaModel(country="SLE",k=10) -> po
plot(po)
coef(po)

stopifnot(
  all.equal(
    coef(po),
    partrans(po,coef(po,transform=TRUE),dir="from")
  ),
  all.equal(
    coef(po,transform=TRUE),
    partrans(po,coef(po),dir="to")
  )
)

pfilter(
  po,
  Np=1000,
  filter.mean=TRUE,
  pred.mean=TRUE,
  pred.var=TRUE,
  filter.traj=TRUE,
  save.states='filter'
) -> pf

logLik(pf)
stopifnot(
  abs(logLik(pf)+100)<0.5
)

plot(pf,yax.flip=TRUE)

simulate(pf) -> sm

plot(cases~day,data=as.data.frame(sm),type="l")
lines(deaths~day,data=as.data.frame(sm),type="l",col="red")

trajectory(po) -> tj
plot(tj,var=c("cases","deaths","I"))

dev.off()
