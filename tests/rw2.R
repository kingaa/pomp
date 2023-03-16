library(pomp)
suppressPackageStartupMessages({
  library(tidyr)
  library(dplyr)
  library(ggplot2)
})

set.seed(1438408329L)

png(filename="rw2-%02d.png",res=100)

rw2(x1_0=1,x2_0=-1) -> po
plot(po)

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
  save.states=TRUE
) -> pf

plot(pf,yax.flip=TRUE)

forecast(pf,format="d") -> fc
simulate(pf) -> sm

emeasure(pf) |> melt() -> ef
vmeasure(pf) |> melt() -> vf
vf |> select(-time,-.id) |> distinct()

bind_rows(
  sim=sm |>
    as.data.frame() |>
    pivot_longer(c(x1,x2,y1,y2)),
  forecast=fc,
  filter=ef |> select(-.id),
  prediction=pred_mean(pf,format="d"),
  filter=filter_mean(pf,format="d"),
  .id="type"
) |>
  ggplot(aes(x=time,y=value,color=factor(type)))+
  geom_line()+
  labs(color="")+
  facet_wrap(~name,scales="free_y")+
  theme_bw()+theme(legend.position="top")

enkf(po,Np=1000) -> kf
plot(kf,yax.flip=TRUE)

eakf(po,Np=1000) -> kf2
plot(kf2,yax.flip=TRUE)

Kf <- kalmanFilter(
  po,
  A=diag(2),
  Q=diag(coef(po,c("s1","s2"))^2),
  C=diag(2),
  R=coef(po,"tau")^2*diag(2)
)

print(c(logLik(pf),logLik(kf),logLik(kf2),Kf$logLik))

stopifnot(
  abs(logLik(pf)+451.6)<0.05,
  abs(logLik(kf)+450.4)<0.05,
  abs(logLik(kf2)+449.1)<0.05,
  abs(Kf$logLik+449.6)<0.05
)

d <- dprocess(sm,log=TRUE)
print(sum(d))
stopifnot(
  abs(sum(d)+381.4) < 0.5
)

dev.off()
