options(digits=3)
png(filename="simulate-%02d.png",res=100)

library(pomp)
library(magrittr)
library(dplyr)

set.seed(1041414791L)

pompExample(ou2)
ou2 %>% simulate(times=0:20,t0=-4,seed=298831503) %>% plot()

try(simulate(rprocess=onestep.sim(Csnippet("z = runif(0,1);")),
  rmeasure=Csnippet("w = rnorm(z,1);"),
  initializer=Csnippet("z = 0;"),
  statenames="z",obsnames="w"))
try(simulate(times=1:10,
  rprocess=onestep.sim(Csnippet("z = runif(0,1);")),
  rmeasure=Csnippet("w = rnorm(z,1);"),
  initializer=Csnippet("z = 0;"),
  statenames="z",obsnames="w"))
try(simulate(times=c(1:5,NA,7:10),t0=0,
  rprocess=onestep.sim(Csnippet("z = runif(0,1);")),
  rmeasure=Csnippet("w = rnorm(z,1);"),
  initializer=Csnippet("z = 0;"),
  statenames="z",obsnames="w"))
try(simulate(times=NULL,t0=0,
  rprocess=onestep.sim(Csnippet("z = runif(0,1);")),
  rmeasure=Csnippet("w = rnorm(z,1);"),
  initializer=Csnippet("z = 0;"),
  statenames="z",obsnames="w"))
try(simulate(times=1:10,t0=NA,
  rprocess=onestep.sim(Csnippet("z = runif(0,1);")),
  rmeasure=Csnippet("w = rnorm(z,1);"),
  initializer=Csnippet("z = 0;"),
  statenames="z",obsnames="w"))
try(simulate(times=1:10,t0=NULL,
  rprocess=onestep.sim(Csnippet("z = runif(0,1);")),
  rmeasure=Csnippet("w = rnorm(z,1);"),
  initializer=Csnippet("z = 0;"),
  statenames="z",obsnames="w"))

simulate(times=1:100,t0=0,seed=450738202,
  rprocess=onestep.sim(Csnippet("z = runif(z-0.5,z+0.5);")),
  rmeasure=Csnippet("w = rnorm(z,1);"),
  initializer=Csnippet("z = 0;"),
  statenames="z",obsnames="w") %>% plot()

try(simulate(times=1:100,t0=0,
  rprocess=onestep.sim(Csnippet("z = runif(z-0.5,z+0.5);")),
  rmeasure=function(t,x,params,...) c(w=rnorm(n=t,x["z"],1)),
  initializer=Csnippet("z = 2;"),
  statenames="z"))

simulate(times=1:100,t0=0,seed=993523767,
  rprocess=onestep.sim(Csnippet("z = runif(z-0.5,z+0.5);")),
  rmeasure=function(t,x,params,...) c(w=rnorm(1,x["z"],1)),
  initializer=Csnippet("z = 0;"),
  statenames="z") -> po
po %>% plot()

simulate(times=1:100,t0=0,seed=378047885,
  rprocess=onestep.sim(function(x,t,params,delta.t,...)
    c(z=runif(n=1,x["z"]-0.5,x["z"]+0.5))),
  rmeasure=function(t,x,params,...) c(w=rnorm(1,x["z"],1)),
  initializer=function(params,t0,...)c(z=0)) %>% plot()

rm(.Random.seed)
simulate(po,seed=406214171) %>% plot()
simulate(po,params=as.list(coef(po)))

data.frame(u=1:10,v=runif(10)) %>%
  pomp(times="u",t0=0) %>%
  simulate(rprocess=onestep.sim(Csnippet("w = runif(0,1);")),
    rmeasure=function(t,x,params,...){
      p <- x["w"]+c(-0.5,0.5)
      c(y=runif(n=1,p[1],p[2]))
    },
    initializer=Csnippet("w=0;"),
    statenames="w") %>%
  obs() %>%
  rownames()

data.frame(u=1:10,v=runif(10)) %>%
  pomp(times="u",t0=0) %>%
  simulate(rprocess=onestep.sim(Csnippet("w = runif(0,1);")),
    rmeasure=Csnippet("y=runif(w-0.5,w+0.5);"),
    initializer=Csnippet("w=0;"),
    statenames="w",obsnames="y") %>%
  obs() %>%
  rownames()

try(simulate(ou2,nsim=-3))
try(simulate(ou2,nsim=NA))
try(simulate(ou2,nsim=NULL))
try(simulate(ou2,nsim="bob"))

ou2 %>% window(end=3) -> po
simulate(po,as.data.frame=TRUE,seed=49569969,nsim=3) %>%
  count(sim)
simulate(po,as.data.frame=TRUE,seed=49569969,nsim=3,include.data=TRUE) %>%
  count(sim)
simulate(po,as.data.frame=TRUE,states=TRUE,seed=49569969)
simulate(po,as.data.frame=TRUE,obs=TRUE,seed=49569969)
simulate(po,as.data.frame=TRUE,obs=FALSE,states=FALSE,seed=49569969)
simulate(po,as.data.frame=TRUE,obs=TRUE,states=TRUE,seed=49569969)
simulate(po,as.data.frame=TRUE,include.data=TRUE,seed=49569969)
simulate(po,as.data.frame=FALSE,include.data=TRUE,seed=49569969)
simulate(po,states=TRUE) %>% rownames()
simulate(po,obs=TRUE) %>% rownames()
simulate(po,obs=TRUE,states=TRUE) %>% names()
simulate(po,nsim=3) %>% show()

data.frame(u=1:10,v=runif(10)) %>%
  pomp(times="u",t0=0) %>%
  simulate(rprocess=onestep.sim(Csnippet("w = runif(0,1);")),
    rmeasure=Csnippet("y=runif(w-0.5,w+0.5);"),
    initializer=Csnippet("w=0;"),
    statenames="w",obsnames="y",as.data.frame=TRUE,include.data=TRUE) -> dat
dat %>% names()
dat %>% dim()

dev.off()
