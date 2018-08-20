options(digits=3)
png(filename="simulate_pomp-%02d.png",res=100)

library(pomp)
library(magrittr)
library(dplyr)

set.seed(1041414791L)

pompExample(ou2)
ou2 -> po

po %>% simulate(times=0:20,t0=-4,seed=298831503) %>% plot()

rm(.Random.seed)
po %>%
  simulate(params=as.list(coef(po)),seed=406214171) %>%
  plot(main="test",yax.flip=TRUE)

set.seed(1041414791L)

data.frame(u=1:10,v=runif(10)) %>%
  pomp(times="u",t0=0) %>%
  simulate(rprocess=onestep.sim(Csnippet("w = runif(0,1);")),
    rmeasure=function(t,x,params,...){
      p <- x["w"]+c(-0.5,0.5)
      c(y=runif(n=1,p[1],p[2]))
    },
    rinit=Csnippet("w=0;"),
    statenames="w") %>%
  obs() %>%
  rownames()

data.frame(u=1:10,v=runif(10)) %>%
  pomp(times="u",t0=0) %>%
  simulate(rprocess=onestep.sim(Csnippet("w = runif(0,1);")),
    rmeasure=Csnippet("y=runif(w-0.5,w+0.5);"),
    rinit=Csnippet("w=0;"),
    statenames="w",obsnames="y") %>%
  obs() %>%
  rownames()

try(simulate(ou2,nsim=-3))
try(simulate(ou2,nsim=NA))
try(simulate(ou2,nsim=NULL))
try(simulate(ou2,nsim="bob"))

ou2 %>% window(end=3) -> po
simulate(po,as.data.frame=TRUE,seed=49569969,nsim=3) %>%
  count(sim) %>% as.data.frame()
simulate(po,as.data.frame=TRUE,seed=49569969,nsim=3,include.data=TRUE) %>%
  count(sim) %>% as.data.frame()
simulate(po,as.data.frame=TRUE,.states=TRUE,seed=49569969)
simulate(po,as.data.frame=TRUE,.obs=TRUE,seed=49569969)
simulate(po,as.data.frame=TRUE,.obs=FALSE,.states=FALSE,seed=49569969)
simulate(po,as.data.frame=TRUE,.obs=TRUE,.states=TRUE,seed=49569969)
simulate(po,as.data.frame=TRUE,include.data=TRUE,seed=49569969)
simulate(po,as.data.frame=FALSE,include.data=TRUE,seed=49569969)
simulate(po,.states=TRUE) %>% rownames()
simulate(po,.obs=TRUE) %>% rownames()
simulate(po,.obs=TRUE,.states=TRUE) %>% names()
simulate(po,nsim=3) %>% show()

data.frame(u=1:10,v=runif(10)) %>%
  pomp(times="u",t0=0) %>%
  simulate(rprocess=onestep.sim(Csnippet("w = runif(0,1);")),
    rmeasure=Csnippet("y=runif(w-0.5,w+0.5);"),
    rinit=Csnippet("w=0;"),
    statenames="w",obsnames="y",as.data.frame=TRUE,include.data=TRUE) -> dat
dat %>% names()
dat %>% dim()

pompExample(rw2)
simulate(rw2,zeronames="x2") %>% plot()

dev.off()
