png(filename="skeleton-%02d.png",res=100)

library(pomp)
library(magrittr)
library(ggplot2)
library(reshape2)
library(plyr)

pompExample(ricker)

ricker <- simulate(ricker,times=1:500,seed=366829807L)
x <- states(ricker)
p <- parmat(coef(ricker),3)
p["r",] <- exp(c(1,2,4))
f <- skeleton(ricker,x=x,params=p,t=time(ricker))
f %>% melt() %>%
  subset(variable=="N",select=-variable) %>%
  dcast(time~rep) %>%
  join(x %>% melt(value.name="x") %>% subset(variable=="N",select=-variable),
       by="time") %>%
  melt(id.vars=c("time","x")) %>%
  mutate(log.r=mapvalues(variable,from=c(1,2,3),to=log(p["r",]))) %>%
  ggplot(aes(x=x,y=value,color=factor(log.r)))+
  geom_line()+
  labs(y=expression(N[t+1]),x=expression(N[t]),color=expression(log(r)))+
  theme_classic()

try(skeleton(x=x,times=time(ricker),params=p))
try(skeleton("ricker",x=x,times=time(ricker),params=p))

pompExample(sir)
p <- parmat(coef(sir),nrep=3)
p["beta2",2:3] <- exp(c(3,5))
trajectory(sir,params=p,times=seq(0,1,length=1000)) -> tj
skeleton(sir,x=tj,params=p,t=seq(0,1,length=1000)) -> f
tj %>% apply(c(1,2),diff) %>% melt(value.name="diff") -> dtj
f %>% melt(value.name="deriv") -> f
join(f,dtj,by=c("time","variable","rep")) %>%
  subset(variable %in% c("S","I","R")) %>%
  ggplot(aes(x=deriv,y=diff/0.001,color=factor(rep)))+
  geom_point()+
  geom_abline(intercept=0,slope=1,color='black')+
  facet_grid(rep~variable,labeller=labeller(rep=label_both))+
  guides(color=FALSE)+
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

dev.off()
