options(digits=3)
png(filename="skeleton-%02d.png",res=100)

library(pomp)
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

ricker() -> ricker

ricker <- simulate(ricker,times=1:500,seed=366829807L)
x <- states(ricker)
p <- parmat(coef(ricker),3)
p["r",] <- exp(c(1,2,4))
f <- skeleton(ricker,x=x,params=p,t=time(ricker))
f %>%
  melt() %>%
  filter(variable=="N") %>%
  select(-variable) %>%
  pivot_wider(names_from=.id) %>%
  left_join(
    x %>% melt() %>% filter(variable=="N") %>%
      select(-variable) %>% rename(x=value),
    by="time"
  ) %>%
  pivot_longer(cols=-c(time,x)) %>%
  mutate(
    name=as.integer(name),
    log.r=log(p["r",name])
  ) %>% 
  ggplot(aes(x=x,y=value,color=factor(signif(log.r,3))))+
  geom_line()+
  labs(y=expression(N[t+1]),x=expression(N[t]),color=expression(log(r)))+
  theme_classic()

skeleton(ricker,params=parmat(coef(ricker),3)) -> dx
stopifnot(dim(dx)==c(2,3,500))

try(skeleton(x=x,times=time(ricker),params=p))
try(skeleton("ricker",x=x,times=time(ricker),params=p))

sir() -> sir
p <- parmat(coef(sir),nrep=3)
p["beta2",2:3] <- exp(c(3,5))
trajectory(sir,params=p,times=seq(0,1,length=1000),format="a") -> tj
skeleton(sir,x=tj,params=p,t=seq(0,1,length=1000)) -> f
tj %>% apply(c(1,2),diff) %>% melt() %>% rename(diff=value) %>%
  mutate(.id=as.integer(.id)) -> dtj
f %>% melt() %>% rename(deriv=value) -> f
right_join(
  f,dtj,
  by=c("time","variable",".id")
) %>%
  filter(variable %in% c("S","I","R")) %>%
  mutate(variable=factor(variable,levels=c("S","I","R"))) %>%
  ggplot(aes(x=deriv,y=diff/0.001,color=factor(.id)))+
  geom_point()+
  geom_abline(intercept=0,slope=1,color='black')+
  facet_grid(.id~variable,labeller=labeller(.id=label_both))+
  guides(color="none")+
  labs(x="derivative",y="finite difference")+
  theme_bw()

try(ricker %>% pomp(skeleton=map(function(...)c(5))) %>%
    skeleton(x=x,times=time(ricker),params=coef(ricker)))
try(ricker %>% pomp(skeleton=map(function(...)c(5,3))) %>%
    skeleton(x=x,times=time(ricker),params=coef(ricker)))
ricker %>% skeleton(x=x,times=time(ricker),params=parmat(coef(ricker),2)) -> xx
try(ricker %>% skeleton(x=xx,times=time(ricker),params=parmat(coef(ricker),3)))
try(ricker %>% skeleton(x=xx,times=time(ricker)[1:5],params=parmat(coef(ricker),2)))

stopifnot(
  ricker %>% pomp(skeleton=NULL) %>%
    skeleton(x=x,times=time(ricker),params=coef(ricker)) %>%
    is.na()
)

try(ricker %>% pomp(skeleton=map(function(...)1,delta.t=-5)))
try(ricker %>% pomp(skeleton=map(function(...)1,delta.t=NA)))
try(ricker %>% pomp(skeleton=map(function(...)1,delta.t=c(1,2,3))))
try(ricker %>% pomp(skeleton=map(function(...)1,delta.t=NULL)))
try(ricker %>% pomp(skeleton=map(function(...)1,delta.t=Inf)))

dev.off()
