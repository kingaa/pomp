options(digits=3)
png(filename="lookup-%02d.png",res=100)

library(pomp)

ct <- covariate_table(x=20:30,y=10:0,times=seq(0,10))
lookup(ct,t=c(1,2.3,4,7))
ct <- covariate_table(x=20:30,y=10:0,times=seq(0,10),order="constant")
lookup(ct,t=c(1,2.3,4,7))
lookup(ct,t=6.1)

try(covariate_table(x=20:30,y=10:0,times=seq(10,0)))
try(covariate_table(x=20:23,y=10:0,times=seq(0,10)))
try(covariate_table(x=20:30,y=10:0,times=seq(0,3)))
covariate_table()
try(covariate_table(times=1:10))
try(covariate_table(a=1:10,times=1))
try(covariate_table(a=1:10,times="a"))
try(covariate_table(a=1:10))
covariate_table(a=1:10,a=10:1,times=1:10)
try(covariate_table(a=1:10,a=10:1,times="a"))
try(covariate_table(data.frame(a=1:10,a=10:1,check.names=FALSE),b=1:10,times="b"))
try(covariate_table(data.frame(a=1:10,a=10:1),b=1:10,times="b"))
try(covariate_table(a=1:10,b=10:1,times="b"))
try(covariate_table(a=1:10,b=10:1,times="c"))
try(covariate_table(a=1:10,b=10:1,times=NA))

library(magrittr)
library(tidyr)
library(ggplot2)

covariate_table(
  bspline.basis(times,nbasis=8,degree=3,deriv=0,names="f0"),
  bspline.basis(times,nbasis=8,degree=3,deriv=1,names="f1"),
  times=seq(0,10,by=0.1)
) %>%
  lookup(t=seq(0,10,by=0.01)) %>%
  gather(variable,value,-t) %>%
  separate(variable,into=c("variable","n")) %>%
  ggplot(aes(x=t,y=value,color=factor(n)))+
  labs(color="element",y="",x="")+
  geom_line()+
  facet_grid(variable~.,
    labeller=labeller(variable=c(f1="derivative",f0="function")))+
  theme_bw()

dev.off()
