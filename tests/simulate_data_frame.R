options(digits=3)

library(pomp)
library(magrittr)
library(dplyr)

set.seed(1041414791L)

data.frame(t=0:10,y=runif(11)) %>%
  simulate(times="t",t0=0,
    rprocess=discrete.time.sim(
      step.fun=function(x,t,params,delta.t,...) {
        setNames(c(params["r"]*x*exp(-x)),"x")
      }),
    rmeasure=function(x,t,params,...) {
      c(y=rpois(n=1,100*x))
    },
    params=c(r=15,x.0=0.1)
  )

try(data.frame(t=0:10,y=runif(11)) %>%
  simulate(times=1:50))
try(data.frame(t=0:10,y=runif(11)) %>%
    simulate(times="t"))
try(data.frame(t=0:10,y=runif(11)) %>%
    simulate(times="t",t0=0))
try(data.frame(t=0:10,y=runif(11)) %>%
    simulate(times="t",t0=0))
try(data.frame(t=0:10,y=runif(11)) %>%
    simulate(times="t",t0=0,params=c(x.0=1)))
