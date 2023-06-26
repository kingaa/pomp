options(digits=3)

library(pomp)

gompertz() -> po

pomp(po,partrans=NULL,bob=3,
  covar=covariate_table(a=0:20,b=0:20,times="a")) -> po1
spy(po1)

po1 |>
  pomp(partrans=parameter_trans(log="r"),params=NULL,
    rinit=function(params,t0,...)params,
    paramnames="r",compile=FALSE,cfile="nancy") -> po2
spy(po2)

sir() -> sir
spy(sir)

rw2() -> rw2
spy(rw2)

sir2() -> sir2
spy(sir2)

try(spy())
try(spy(list()))

pomp(
  data=NULL,
  t0=0,times=1:10,
  x0=as.double(1),
  params=c(x_0=1,a=22),
  rinit=Csnippet(r"{x = *get_userdata_double("x0");}"),
  statenames="x",
  compile=FALSE
) |>
  spy()
