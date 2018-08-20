options(digits=3)
png(filename="simulate_missing-%02d.png",res=100)

library(pomp)
library(magrittr)
library(dplyr)

set.seed(1041414791L)

try(simulate(rprocess=onestep.sim(Csnippet("z = runif(0,1);")),
  rmeasure=Csnippet("w = rnorm(z,1);"),
  rinit=Csnippet("z = 0;"),
  statenames="z",obsnames="w"))
try(simulate(times=1:10,
  rprocess=onestep.sim(Csnippet("z = runif(0,1);")),
  rmeasure=Csnippet("w = rnorm(z,1);"),
  rinit=Csnippet("z = 0;"),
  statenames="z",obsnames="w"))
try(simulate(times=c(1:5,NA,7:10),t0=0,
  rprocess=onestep.sim(Csnippet("z = runif(0,1);")),
  rmeasure=Csnippet("w = rnorm(z,1);"),
  rinit=Csnippet("z = 0;"),
  statenames="z",obsnames="w"))
try(simulate(times=NULL,t0=0,
  rprocess=onestep.sim(Csnippet("z = runif(0,1);")),
  rmeasure=Csnippet("w = rnorm(z,1);"),
  rinit=Csnippet("z = 0;"),
  statenames="z",obsnames="w"))
try(simulate(times=1:10,t0=NA,
  rprocess=onestep.sim(Csnippet("z = runif(0,1);")),
  rmeasure=Csnippet("w = rnorm(z,1);"),
  rinit=Csnippet("z = 0;"),
  statenames="z",obsnames="w"))
try(simulate(times=1:10,t0=NULL,
  rprocess=onestep.sim(Csnippet("z = runif(0,1);")),
  rmeasure=Csnippet("w = rnorm(z,1);"),
  rinit=Csnippet("z = 0;"),
  statenames="z",obsnames="w"))

simulate(times=1:100,t0=0,seed=450738202,
  rprocess=onestep.sim(Csnippet("z = runif(z-0.5,z+0.5);")),
  rmeasure=Csnippet("w = rnorm(z,1);"),
  rinit=Csnippet("z = 0;"),
  statenames="z",obsnames="w") %>% plot()

try(simulate(times=1:100,t0=0,
  rmeasure=Csnippet("w = rnorm(z,1);"),
  rinit=Csnippet("z = 0;"),
  statenames="z",obsnames="w"))

simulate(times=1:100,t0=0,
  rprocess=onestep.sim(function(x,t,params,delta.t,...) {
    z <- x["z"]
    c(z=runif(1,z-0.5,z+0.5))
  }),
  rmeasure=function(x,t,params,...) {
    z <- x["z"]
    c(w = rnorm(1,z,1))
  },params=c(z.0=0))

try(simulate(times=1:100,t0=0,
  rprocess=onestep.sim(function(x,t,params,delta.t,...) {
    z <- x["z"]
    c(z=runif(1,z-0.5,z+0.5))
  }),params=c(z.0=0)))

simulate(times=1:100,t0=0,
  rprocess=onestep.sim(function(x,t,params,delta.t,...) {
    z <- x["z"]
    c(z=runif(1,z-0.5,z+0.5))
  }),params=c(z.0=0),.states=TRUE) %>%
  str()

try(simulate(times=1:100,t0=0,
  rprocess=onestep.sim(function(x,t,params,delta.t,...) {
    z <- x["z"]
    c(z=runif(1,z-0.5,z+0.5))
  }),params=c(z.0=0),.obs=TRUE))

simulate(times=1:100,t0=0,
  rprocess=onestep.sim(function(x,t,params,delta.t,...) {
    z <- x["z"]
    c(z=runif(1,z-0.5,z+0.5))
  }),
  rmeasure=function(x,t,params,...) {
    z <- x["z"]
    c(w = rnorm(1,z,1))
  },
  params=c(z.0=0),.obs=TRUE) %>%
  str()

simulate(times=1:100,t0=0,
  rprocess=onestep.sim(function(x,t,params,delta.t,...) {
    z <- x["z"]
    c(z=runif(1,z-0.5,z+0.5))
  }),
  rmeasure=function(x,t,params,...) {
    z <- x["z"]
    c(w = rnorm(1,z,1))
  },
  params=c(z.0=0),.states=TRUE,.obs=TRUE) %>%
  str()

try(simulate(times=1:100,t0=0,
  rprocess=onestep.sim(Csnippet("z = runif(z-0.5,z+0.5);")),
  rmeasure=function(t,x,params,...) c(w=rnorm(n=t,x["z"],1)),
  rinit=Csnippet("z = 2;"),
  statenames="z"))

simulate(times=1:100,t0=0,seed=993523767,
  rprocess=onestep.sim(Csnippet("z = runif(z-0.5,z+0.5);")),
  rmeasure=function(t,x,params,...) c(w=rnorm(1,x["z"],1)),
  rinit=Csnippet("z = 0;"),
  statenames="z") %>% plot()

simulate(times=1:100,t0=0,seed=378047885,
  rprocess=onestep.sim(function(x,t,params,delta.t,...)
    c(z=runif(n=1,x["z"]-0.5,x["z"]+0.5))),
  rmeasure=function(t,x,params,...) c(w=rnorm(1,x["z"],1)),
  rinit=function(params,t0,...)c(z=0)) %>% plot()

dev.off()
