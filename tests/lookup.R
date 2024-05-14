options(digits=3)
suppressPackageStartupMessages({
  library(tidyr)
  library(ggplot2)
})

png(filename="lookup-%02d.png",res=100)

library(pomp)

ct <- covariate_table(x=20:30,y=10:0,times=seq(0,10))
lookup(ct,t=c(1,2.3,4,7,20))
ct <- covariate_table(x=20:30,y=10:0,times=seq(0,10),order="constant")
lookup(ct,t=c(1,2.3,4,7))
lookup(ct,t=6.1)
plot(y~t,data=lookup(ct,t=seq(0,10.5,by=0.01)),type='l')
lines(seq(0,10),10:0,col="blue",type="s")
plot(x~t,data=lookup(ct,t=seq(0,10.5,by=0.01)),type='l')
lines(seq(0,10),20:30,col="blue",type="s")

ct <- covariate_table(x=20:31,y=12:1,times=c(0:5,5:10),order="constant")
lookup(ct,t=seq(4.9,5.1,by=0.05))

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

covariate_table(
  bspline_basis(times,nbasis=8,degree=3,deriv=0,names="f0"),
  bspline_basis(times,nbasis=8,degree=3,deriv=1,names="f1"),
  times=seq(0,10,by=0.1)
) |>
  lookup(t=seq(0,10,by=0.01)) |>
  pivot_longer(-t) |>
  separate(name,into=c("variable","n")) |>
  ggplot(aes(x=t,y=value,color=factor(n)))+
  labs(color="element",y="",x="")+
  geom_line()+
  facet_grid(variable~.,
    labeller=labeller(variable=c(f1="derivative",f0="function")))+
  theme_bw()

covariate_table(
  x=seq(1,10,by=1),
  y=seq(2,20,by=2),
  times="x"
) -> tab
lookup(tab,c(0,2.5,3.6,10,20))

repair_lookup_table(tab,t=c(seq(0,10),20)) -> tab2
lookup(tab2,c(0,2.5,3.6,10,20))

covariate_table(
  x=seq(1,10,by=1),
  y=seq(1,10,by=1),
  order="const",
  times="x"
) -> tab
lookup(tab,c(0,2.5,3.6,10,20))

repair_lookup_table(tab,t=c(seq(0,10,by=1),20)) -> tab2
lookup(tab2,c(0,2.5,3.6,10,20))

covariate_table(
  x=seq(1,10,by=1),
  y=seq(1,10,by=1),
  order="const",
  times="x"
) -> tab

repair_lookup_table(tab,t=c(seq(0,10,by=1),20),order="lin") -> tab2
lookup(tab2,c(0,2.5,3.6,10,20))

repair_lookup_table(tab,t=c(seq(0,10,by=1),20),order="const") -> tab2
lookup(tab2,c(0,2.5,3.6,10,20))

dev.off()
