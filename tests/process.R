options(digits=3)

library(pomp)

ou2() |> window(end=10) -> po

set.seed(293982095)

po |>
  simulate(format="arrays",nsim=3,seed=3434388L) |>
  getElement("states") -> x
po |> time() -> t
coef(po) |> parmat(7) -> p
p["sigma_1",] <- seq(from=1,to=7,by=1)

try(dprocess("ou2",x=x,times=t,params=p))
try(dprocess(x=x,times=t,params=p))
try(po |> dprocess(x=x,times=t,params=p,log=TRUE))
try(po |> dprocess(x=x,times=t,params=p[,1:2],log=TRUE))
try(po |> dprocess(x=x,times=t[1:5],params=p,log=TRUE))
po |> dprocess(x=x,times=t,params=p[,1:3],log=TRUE) -> d1
po |> dprocess(x=x,times=t,params=p[,2],log=TRUE) -> d2
stopifnot(d1[2,]==d2[2,])
try(po |> dprocess(x=x[,,2],times=t[2],params=p[,1:3],log=FALSE))
try(po |> dprocess(x=x[,,2:5],times=t[2:5],params=p[,1:2],log=FALSE))
po |> dprocess(x=x[,,2:5],times=t[2:5],params=p[,1:3],log=TRUE) |>
  apply(1,sum)
try(po |> dprocess(x=x[,1:2,2:5],times=t[2:5],params=p[,1:3],log=TRUE) |>
      apply(1,sum))
po |> dprocess(x=x[,1,2:5],times=t[2:5],params=p[,1:3],log=TRUE) |>
  apply(1,sum)
po |> pomp(dprocess=NULL) |>
  dprocess(x=x[,1,2:5],times=t[2:5],params=p[,1:3],log=TRUE) |> is.na() |>
  stopifnot()

po |> rinit(params=coef(po)) -> x0
freeze(po |>
         rprocess(params=coef(po),x0=parmat(x0,3),t0=timezero(po),times=time(po),
           offset=1),
  seed=3434388L) -> x1

stopifnot(max(abs(x-x1))==0)

po |> rinit(nsim=6) -> x0
try(rprocess("ou2",x0=x0,t0=t[1],times=t,params=p))
try(rprocess(x0=x0,t0=t[1],times=t,params=p))
freeze(po |> rprocess(times=t,params=p),seed=995484) -> x1
try(po |> rprocess(x0=x0,params=p))
try(po |> rprocess(x0=x0,t0=t[1],params=p))
freeze(po |> rprocess(x0=x0,times=t),seed=995484) -> x2
try(po |> rprocess(x0=x0,times=t,params=p))
try(po |> rprocess(x0=x0,t0=t[1],times=t,params=p))
po |> rprocess(x0=x0,t0=t[1],times=t,params=p[,1:3]) -> x
stopifnot(
  dim(x)==c(2,6,10),
  names(dimnames(x))==c("variable",".id","time")
)
po |> rprocess(x0=x0[,2],t0=t[1],times=t,params=p[,1:3]) -> x
stopifnot(
  dim(x)==c(2,3,10),
  names(dimnames(x))==c("variable",".id","time")
)

try(po |> rprocess(x0=x0,t0=t[1],times=numeric(0),params=p))
try(po |> rprocess(x0=x0[,2],t0=t[2],times=t[1],params=p[,1:3]))
try(po |> rprocess(x0=x0[,2:4],t0=t[2],times=t[2:5],params=p[,1:2]))
po |> rprocess(x0=x0[,2:4],t0=t[2],times=t[2:5],params=p[,1:3]) |>
  apply(1,sum)

simulate(
  times=seq(0,10), t0=0,
  params=c(s=3,x_0=0,tau=1),
  covar=covariate_table(z=c(1,1),times=c(0,10)),
  rprocess = onestep(
    function (t, x, s, delta.t, ...) {
      c(x=rnorm(n=1,mean=x,sd=s*sqrt(delta.t)))
    }
  ),
  dprocess = function (x_1, x_2, s, t_1, t_2, ...) {
    delta.t <- t_2-t_1
    dnorm(x=x_2,mean=x_1,sd=s*sqrt(delta.t),log=TRUE)
  },
  seed=3434388L
) -> rw

dprocess(rw,x=states(rw),params=coef(rw),times=time(rw),log=TRUE) -> d
stopifnot(
  round(sum(d),1)==-23.2,
  dim(d)==c(1,10),
  names(dimnames(d))==c(".id","time")
)
