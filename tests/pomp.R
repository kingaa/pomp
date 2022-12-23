options(digits=3)

set.seed(758723694)

library(pomp)

try(pomp())
try(pomp("bob"))
try(pomp(times=3))
try(pomp(NULL))
try(data.frame(a=1:10,a=1:10,check.names=FALSE) |> pomp(t0=4))
try(data.frame(a=1:10,a=1:10,check.names=FALSE) |> pomp(times=1:10,t0=4))
try(data.frame(a=1:10,a=1:10,check.names=FALSE) |> pomp(times="a",t0=4))
try(data.frame(a=1:10,b=1:10) |> pomp(times=1:10,t0=4))
try(data.frame(a=1:10,b=1:10) |> pomp(times="b"))
try(data.frame(a=10:1,b=1:10) |> pomp(times="a",t0=0))
data.frame(a=1:10,b=1:10) |> pomp(times=1:10,t0=0)
try(data.frame(a=1:10,b=1:10) |> pomp(times=1))
try(data.frame(a=1:10,b=1:10) |> pomp(times="a"))
try(data.frame(a=1:10,b=1:10) |> pomp(times="c",t0=0))
try(data.frame(a=1:10,b=1:10) |> pomp(times=NA,t0=0))
try(data.frame(a=1:10,b=1:10) |> pomp(times=NULL,t0=0))
try(data.frame(a=1:10,b=1:10) |> pomp(times=1,t0=11))
try(data.frame(a=1:10,b=1:10) |> pomp(times="a",t0=11))
try(data.frame(a=1:10,b=1:10) |> pomp(times="b",t0=NULL))
try(data.frame(a=1:10,b=1:10) |> pomp(times="a",t0=NA))
stopifnot(data.frame(a=1:10,b=1:10) |>
    pomp(covar=covariate_table(c=0:10,d=0:10,times="c"),
      covarnames="d",times="a",t0=0,bob=3) |> class() == "pomp")
try(data.frame(a=1:10,b=1:10) |>
    pomp(covar=covariate_table(c=1:10,d=1:10,d=1:10,times="c"),
      times="a",t0=0))
stopifnot(data.frame(a=1:10,b=1:10) |>
    pomp(times="a",t0=0) |> class() == "pomp")
try(NULL |> pomp(t0=4))
try(NULL |> pomp(times="a",t0=0))
try(NULL |> pomp(times=1:10,t0=3))
try(NULL |> pomp(times=1:10,t0=1,rinit=3))
stopifnot(NULL |> pomp(times=1:10,t0=1) |> class() == "pomp")

gompertz() -> po
stopifnot({
  po |> pomp(rprocess=NULL) -> po1
  rprocess(po1,x0=rinit(po1),t0=timezero(po1),
    times=time(po1),params=coef(po1)) -> x
  sum(is.na(x))==100
})

po |> pomp(rprocess=NULL) |> slot("rprocess")
po |> pomp(skeleton=NULL) |> slot("skeleton")
po |> pomp(partrans=NULL) |> slot("partrans")

stopifnot({
  po |> pomp(partrans=NULL) |> coef(transform=TRUE) -> theta1
  coef(po) -> theta2
  theta1==theta2
},
  po |> pomp(times=1:5) |> class() == "pomp")

stopifnot(po |>
    pomp(rprocess=onestep(function(x,t,params,delta.t,...)x),
      skeleton=map(function(x,t,params,...)x),
      rmeasure=function(...)3,
      dmeasure=function(log,...)1,
      covar=covariate_table(a=1:20,b=1:20,times="a")) |> class() == "pomp")

try(po |> pomp(times=3:1))

try(po |> pomp(rinit=Csnippet("X=3;")))
stopifnot(
  po |> pomp(rinit=Csnippet("X=3;"),statenames=c("X","Z")) |>
    class() == "pomp")

try(po |> pomp(rprocess="bob"))
try(po |> pomp(skeleton="bob"))
try(po |> pomp(partrans="bob"))
try(po |> pomp(params=c(1,2,3)))
try(po |> pomp(params=c(a=1,b=2,3)))

sir() |> window(end=0.12) -> po2
po2 |> simulate(seed=4358686) |> as.data.frame()
pomp(po2,covar=NULL)@covar
try(po2 |> pomp(covar="bob"))
try(po2 |> pomp(rmeasure=function(x)x))
try(pomp(data=NULL,times=1:10,t0=0,rmeasure=Csnippet("")))

try(po2 |> pomp(rmeasure=Csnippet("reports=3;"),cfile="sir_source"))
