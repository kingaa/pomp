options(digits=3)
png(filename="simulate_data_frame-%02d.png",res=100)

library(pomp)
library(magrittr)
library(dplyr)

set.seed(1041414791L)

data.frame(t=0:10,y=runif(11)) -> dat

dat %>%
  simulate(times="t",t0=0,
    rprocess=discrete.time.sim(
      step.fun=function(x,t,params,delta.t,...) {
        setNames(c(params["r"]*x*exp(-x)),"x")
      }),
    rmeasure=function(x,t,params,...) {
      c(y=rpois(n=1,100*x))
    },
    params=c(r=15,x.0=0.1)
  ) %>% plot()

try(dat %>% simulate(times=1:50))
try(dat %>% simulate(times="t"))
try(dat %>% simulate(times="t",t0=0))
try(dat %>% simulate(times=1,t0=0))
try(dat %>% simulate(times="t",t0=NA))
try(dat %>% simulate(times=NA,t0=0))
try(dat %>% simulate(times="t",t0=0,params=c(x.0=1)))

dev.off()
