options(digits=3)
png(filename="simulate-%02d.png",res=100)

library(pomp)
library(dplyr)

set.seed(1041414791L)

ou2() -> ou2
ou2 |> simulate(times=0:20,t0=-4,seed=298831503) |> plot()

try(simulate(rprocess=onestep(Csnippet("z = runif(0,1);")),
  rmeasure=Csnippet("w = rnorm(z,1);"),
  rinit=Csnippet("z = 0;"),
  statenames="z",obsnames="w"))
try(simulate(times=1:10,
  rprocess=onestep(Csnippet("z = runif(0,1);")),
  rmeasure=Csnippet("w = rnorm(z,1);"),
  rinit=Csnippet("z = 0;"),
  statenames="z",obsnames="w"))
try(simulate(times=c(1:5,NA,7:10),t0=0,
  rprocess=onestep(Csnippet("z = runif(0,1);")),
  rmeasure=Csnippet("w = rnorm(z,1);"),
  rinit=Csnippet("z = 0;"),
  statenames="z",obsnames="w"))
try(simulate(times=NULL,t0=0,
  rprocess=onestep(Csnippet("z = runif(0,1);")),
  rmeasure=Csnippet("w = rnorm(z,1);"),
  rinit=Csnippet("z = 0;"),
  statenames="z",obsnames="w"))
try(simulate(times=1:10,t0=NA,
  rprocess=onestep(Csnippet("z = runif(0,1);")),
  rmeasure=Csnippet("w = rnorm(z,1);"),
  rinit=Csnippet("z = 0;"),
  statenames="z",obsnames="w"))
try(simulate(times=1:10,t0=NULL,
  rprocess=onestep(Csnippet("z = runif(0,1);")),
  rmeasure=Csnippet("w = rnorm(z,1);"),
  rinit=Csnippet("z = 0;"),
  statenames="z",obsnames="w"))

simulate(times=1:100,t0=0,seed=450738202,
  rprocess=onestep(Csnippet("z = runif(z-0.5,z+0.5);")),
  rmeasure=Csnippet("w = rnorm(z,1);"),
  rinit=Csnippet("z = 0;"),
  statenames="z",obsnames="w") |> plot()

try(simulate(times=1:100,t0=0,
  rmeasure=Csnippet("w = rnorm(z,1);"),
  rinit=Csnippet("z = 0;"),
  statenames="z",obsnames="w"))

simulate(times=1:100,t0=0,
  rprocess=onestep(function(z,...) {
    c(z=runif(1,z-0.5,z+0.5))
  }),
  rmeasure=function(z,...) {
    c(w = rnorm(1,z,1))
  },params=c(z.0=0))

stopifnot(
  {
    simulate(times=1:100,t0=0,
      rprocess=onestep(function(z,...) {
        c(z=runif(1,z-0.5,z+0.5))
      }),params=c(z.0=0)) |> obs() -> y
    dim(y)==c(0,100)
  }
)

stopifnot(
  {
    simulate(times=1:100,t0=0,
      rprocess=onestep(function(z,...) {
        c(z=runif(1,z-0.5,z+0.5))
      }),params=c(z.0=0),format="arrays") -> s
    dim(s$states)==c(1,1,100)
    dim(s$obs)==c(0,1,100)
  }
)

try(simulate(times=1:100,t0=0,
  rprocess=euler(function(w,z,...,delta.t) c(w,z),delta.t=0.1),
  params=c(w_0=1,z_0=3)))

try(simulate(times=1:100,t0=0,
  rprocess=euler(function(w,z,...,delta.t) c(w=w,z),delta.t=0.1),
  params=c(w_0=1,z_0=3)))

try(simulate(times=1:100,t0=0,
  rprocess=onestep(Csnippet("z = runif(z-0.5,z+0.5);")),
  rmeasure=function(t,z,...) c(w=rnorm(n=t,z,1)),
  rinit=Csnippet("z = 2;"),
  statenames="z"))

simulate(times=1:100,t0=0,seed=993523767,
  rprocess=onestep(Csnippet("z = runif(z-0.5,z+0.5);")),
  rmeasure=function(z,...) c(w=rnorm(1,z,1)),
  rinit=Csnippet("z = 0;"),
  statenames="z") -> po
po |> plot()

simulate(times=1:100,t0=0,seed=378047885,
  rprocess=onestep(function(z,...)
    c(z=runif(n=1,z-0.5,z+0.5))),
  rmeasure=function(z,...) c(w=rnorm(1,z,1)),
  rinit=function(params,t0,...)c(z=0)) |> plot()

rm(.Random.seed)
po |>
  simulate(params=as.list(coef(po)),seed=406214171) |>
  plot(variables=rep(c("z","w"),10),main="test",yax.flip=TRUE)

set.seed(1041414791L)

data.frame(u=1:10,v=runif(10)) |>
  pomp(times="u",t0=0) |>
  simulate(rprocess=onestep(Csnippet("w = runif(0,1);")),
    rmeasure=function(w,...){
      p <- w+c(-0.5,0.5)
      c(y=runif(n=1,p[1],p[2]))
    },
    rinit=Csnippet("w=0;"),
    statenames="w") |>
  obs() |>
  rownames()

data.frame(u=1:10,v=runif(10)) |>
  pomp(times="u",t0=0) |>
  simulate(rprocess=onestep(Csnippet("w = runif(0,1);")),
    rmeasure=Csnippet("y=runif(w-0.5,w+0.5);"),
    rinit=Csnippet("w=0;"),
    statenames="w",obsnames="y") |>
  obs() |>
  rownames()

try(simulate(ou2,nsim=-3))
try(simulate(ou2,nsim=NA))
try(simulate(ou2,nsim=NULL))
try(simulate(ou2,nsim="bob"))

ou2 |> window(end=3) -> po
simulate(po,format="data.frame",seed=49569969,nsim=3) |>
  count(.id) |> as.data.frame()
simulate(po,format="data.frame",seed=49569969,nsim=3,include.data=TRUE) |>
  count(.id) |> as.data.frame()
simulate(po,format="data.frame",seed=49569969)
simulate(po,format="data.frame",include.data=TRUE,seed=49569969)
simulate(po,format="arrays") -> s
s$states |> rownames()
s$obs |> rownames()
simulate(po,nsim=3) |> show()

data.frame(u=1:10,v=runif(10)) |>
  pomp(times="u",t0=0) |>
  simulate(rprocess=onestep(Csnippet("w = runif(0,1);")),
    rmeasure=Csnippet("y=runif(w-0.5,w+0.5);"),
    rinit=Csnippet("w=0;"),
    statenames="w",obsnames="y",format="data.frame",include.data=TRUE) -> dat
dat |> names()
dat |> dim()

rw2() -> rw2
simulate(rw2,accumvars="x2") |> plot()
try(simulate(rw2,params=c(a="bob",b="nancy")))

simulate(rw2,times=c(1:5,5,5),format="d",seed=49569969) -> x
stopifnot(with(x,c(x1[6:7]==x1[5],x2[6:7]==x2[5],time[6:7]==time[5],
  y1[6:7]!=y1[5],y2[6:7]!=y2[5])))

simulate(rw2,rinit=function(...) c(x1=0,x2=0)) -> x

library(tidyr)
dacca() |>
  simulate(nsim=2,include.data=TRUE,format="data.frame") -> dat
dat |>
  select(time,.id,seas_1,seas_2) |>
  arrange(time,.id) |>
  head(6)
dat |>
  pivot_longer(-.id) |>
  group_by(.id,name) |>
  summarize(n=sum(is.na(value)),.groups="drop") |>
  filter(n>0) -> d1
stopifnot(
  `na in sim`=all(d1$.id=="data"),
  full=all(d1$n==600)
)

dev.off()
