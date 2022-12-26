library(pomp)
suppressPackageStartupMessages({
  library(tidyr)
  library(dplyr)
  library(ggplot2)
})

set.seed(1438408329L)

png(filename="ou2-%02d.png",res=100)

ou2() -> po
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

emeasure(pf) %>% melt() -> ef
vmeasure(pf) %>% melt() -> vf
vf %>% select(-time,-.id) %>% distinct()

bind_rows(
  sim=pf %>%
    as.data.frame() %>%
    pivot_longer(c(x1,x2,y1,y2),names_to="variable"),
  forecast=fc,
  filter=ef %>% select(-.id),
  prediction=pred_mean(pf,format="d"),
  filter=filter_mean(pf,format="d"),
  .id="type"
) %>%
  ggplot(aes(x=time,y=value,color=factor(type)))+
  geom_line()+
  labs(color="")+
  facet_wrap(~variable,scales="free_y")+
  theme_bw()+theme(legend.position="top")

enkf(po,Np=1000) -> kf
plot(kf,yax.flip=TRUE)

eakf(po,Np=1000) -> kf2
plot(kf2,yax.flip=TRUE)

Kf <- kalmanFilter(
  po,
  A=matrix(coef(po,c("alpha_1","alpha_2","alpha_3","alpha_4")),2,2),
  Q={
    q <- matrix(c(coef(po,c("sigma_1","sigma_2")),0,coef(po,"sigma_3")),2,2)
    tcrossprod(q)
  },
  C=diag(2),
  R=diag(coef(po,"tau")^2,2)
)

stopifnot(
  abs(logLik(pf)+478.7)<0.5,
  abs(logLik(kf)+477.1)<0.5,
  abs(logLik(kf2)+476.3)<0.5,
  abs(Kf$logLik+476.6)<0.5
)

trajectory(po) -> tj
plot(tj)

d <- dprocess(sm,log=TRUE)
stopifnot(
  abs(sum(d)+452.5)<0.5
)

dev.off()
