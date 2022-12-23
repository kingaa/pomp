library(pomp)

try(emeasure())
try(emeasure("does this work?"))
try(simulate(t0=0,times=1:10,emeasure=Csnippet("E_B=rho*inc;")))

ou2() -> po
emeasure(po) -> x1
emeasure(po,x=states(po)) -> x2
emeasure(po,x=states(po),params=coef(po)) -> x3
stopifnot(
  dim(x1)==c(2,1,100),
  identical(x1,x2),
  identical(x1,x3)
)
try(emeasure(po,x=states(po),params=coef(po),times=numeric(0)))
try(emeasure(po,x=states(po),params=coef(po),times=c(1,2,3)))
simulate(po,nsim=3,format="arrays") |> getElement("states") -> X
try(emeasure(po,x=X,params=parmat(coef(po),2),times=time(po)))
po |> simulate(emeasure=function(x1, x2, ...) x1+x2) -> po1
try(emeasure(po1,x=states(po1),params=coef(po1),times=time(po1)))
po |> simulate(emeasure=function(x1, x2, t, ...)
  setNames(rep(x1+x2,ceiling(t)),head(letters,ceiling(t)))) -> po1
try(emeasure(po1,x=states(po1),params=coef(po1),times=time(po1)))
po |> simulate(emeasure=NULL) -> po1
e <- emeasure(po1,x=states(po1),params=coef(po1),times=time(po1))
stopifnot(
  dim(e)==c(2,1,100),
  sum(is.na(e))==200
)

sir() |>
  simulate(
    times=(1:10)/52,
    emeasure=function(cases, rho, seas_1, seas_2, seas_3, ...)
      c(reports=cases*rho)
  ) -> po
e <- emeasure(po,x=states(po)[,10],
  params=coef(po),times=time(po)[10])
