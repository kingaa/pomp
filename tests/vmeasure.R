library(pomp)

try(vmeasure())
try(vmeasure("does this work?"))
try(simulate(t0=0,times=1:10,vmeasure=Csnippet("V_B_B=rho*(1-rho)*inc;")))

ou2() -> po
vmeasure(po) -> v1
vmeasure(po,x=states(po)) -> v2
vmeasure(po,x=states(po),params=coef(po)) -> v3
stopifnot(
  dim(v1)==c(2,2,1,100),
  identical(v1,v2),
  identical(v1,v3)
)
try(vmeasure(po,x=states(po),params=coef(po),times=numeric(0)))
try(vmeasure(po,x=states(po),params=coef(po),times=c(1,2,3)))
simulate(po,nsim=3,format="arrays") |> getElement("states") -> X
try(vmeasure(po,x=X,params=parmat(coef(po),2),times=time(po)))
po |>
  simulate(
    vmeasure=function(x1, x2, ...) {
      matrix(c(x1+x2,1,1,x1+x2,1,1),nrow=2)
    }
  ) -> po1
try(vmeasure(po1,x=states(po1),params=coef(po1),times=time(po1)))
po |>
  simulate(
    vmeasure=function(x1, x2, ...) {
      matrix(c(x1+x2,1,1,x1+x2),nrow=2)
    }
  ) -> po1
try(vmeasure(po1,x=states(po1),params=coef(po1),times=time(po1)))
po |>
  simulate(
    vmeasure=function(x1, x2, t, ...) {
      v <- matrix(rep(x1+x2,ceiling(t)^2),nrow=ceiling(t))
      nm <- head(letters,ceiling(t))
      dimnames(v) <- list(nm,nm)
      v
    }
  ) -> po1
try(vmeasure(po1,x=states(po1),params=coef(po1),times=time(po1)))
po |>
  simulate(
    vmeasure=function(x1, x2, t, ...) {
      v <- matrix(c(x1+x2,1,1,x1+x2),nrow=2)
      nm <- c("a","b")
      dimnames(v) <- list(nm,NULL)
      v
    }
  ) -> po1
v <- vmeasure(po1,x=states(po1),params=coef(po1),times=time(po1))
stopifnot(
  dim(v)==c(2,2,1,100),
  v[1,1,,]==v[2,2,,],
  v[1,2,,]==v[2,1,,],
  rownames(v)==c("a","b"),
  colnames(v)==c("a","b"),
  sum(is.na(v))==0
)
po |> simulate(vmeasure=NULL) -> po1
v <- vmeasure(po1,x=states(po1),params=coef(po1),times=time(po1))
stopifnot(
  dim(v)==c(2,2,1,100),
  rownames(v)==c("y1","y2"),
  colnames(v)==c("y1","y2"),
  sum(is.na(v))==400
)

sir() |>
  simulate(
    times=(1:10)/52,
    vmeasure=function(cases, rho, seas_1, seas_2, seas_3, ...) {
      v <- array(dim=c(1,1),dimnames=list("reports","reports"))
      v[,] <- cases*rho*(1-rho)
      v
    }
  ) -> po
v <- vmeasure(po,x=states(po)[,10],params=coef(po),times=time(po)[10])
