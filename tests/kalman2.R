options(digits=3)
png(filename="kalman2-%02d.png",res=100)

library(pomp)
suppressPackageStartupMessages({
  library(tidyr)
})

set.seed(1968726372)

gompertz() -> po

po |>
  window(start=1,end=30) |>
  as.data.frame() |>
  subset(select=-X) -> dat

try(dat |> enkf())

dat |>
  enkf(
    times="time",t0=0,Np=100,
    params=c(r=0.1,K=150),
    rinit=function(K, ...) {
      c(x=K)
    },
    rprocess=discrete_time(
      function (x, r, K, ...) {
        e <- rnorm(n=1,mean=0,sd=0.1)
        c(x=x*exp(r*(1-x/K))+e)
      }
    ),
    emeasure=function(x,...) c(Y=0.01*x),
    vmeasure=function(...) matrix(2,1,1,dimnames=list("Y","Y"))
  ) |> plot()

try(dat |> eakf())

dat |>
  eakf(
    times="time",t0=0,Np=100,
    params=c(r=0.1,K=100),
    rinit=function(K, ...) {
      c(x=rlnorm(n=1,meanlog=log(K),sdlog=0.01))
    },
    rprocess=discrete_time(
      function (x, r, K, ...) {
        c(x=x*exp(r*(1-x/K)))
      }
    ),
    vmeasure=function(R, ...) R,
    emeasure=function(x,...) c(Y=0.01*x),
    R=matrix(0.01,1,1,dimnames=list("Y","Y"))
  ) -> kf

kf |> plot()

kf |> as.data.frame() |> names()
kf |> as.data.frame() |> pivot_longer(-time) |> names()
kf |> forecast() |> melt() |> sapply(class)
kf |> forecast(format="d") |> sapply(class)
kf |> filter_mean() |> melt() |> sapply(class)
kf |> filter_mean(format="d") |> sapply(class)
kf |> pred_mean() |> melt() |> sapply(class)
kf |> pred_mean(format="d") |> sapply(class)
try(kf |> pred_var() |> melt() |> sapply(class))
try(kf |> pred_var(format="d") |> sapply(class))
try(kf |> filter_traj() |> melt() |> sapply(class))
try(kf |> filter_traj(format="d") |> sapply(class))
try(kf |> saved_states() |> melt() |> sapply(class))
try(kf |> saved_states(format="d") |> sapply(class))
try(kf |> forecast(format="l"))

kf |> forecast(format="a") |> melt() -> kdat1a
kf |> filter_mean(format="a") |> melt() -> kdat2a
kf |> pred_mean(format="a") |> melt() -> kdat3a
try(kf |> filter_traj(format="a") |> melt() -> kdat4a)
try(kf |> eff_sample_size(format="n") |> melt() -> kdat5a)
kf |> cond_logLik(format="n") |> melt() -> kdat6a
try(kf |> saved_states(format="l") |> melt() -> kdat7a)
kf |> as.data.frame() -> kdat0
kf |> forecast(format="d") -> kdat1
kf |> filter_mean(format="d") -> kdat2
kf |> pred_mean(format="d") -> kdat3
try(kf |> filter_traj(format="d") -> kdat4)
try(kf |> eff_sample_size(format="d") -> kdat5)
kf |> cond_logLik(format="d") -> kdat6
try(kf |> saved_states(format="d") -> kdat7)
stopifnot(
  all(kdat0$forecast.Y==kdat1$value),
  all(kdat0$filter.mean.x==kdat2$value),
  all(kdat0$pred.mean.x==kdat3$value),
  all(kdat0$cond.logLik==kdat6$value),
  all.equal(kdat1$time,as.numeric(kdat1a$time)),
  all.equal(kdat2$time,as.numeric(kdat2a$time)),
  all.equal(kdat3$time,as.numeric(kdat3a$time)),
  all.equal(kdat1$value,kdat1a$value),
  all.equal(kdat2$value,kdat2a$value),
  all.equal(kdat3$value,kdat3a$value),
  all.equal(kdat6$cond.logLik,kdat6a$value)
)

try(po |> enkf(rprocess=NULL))
try(po |> eakf(rprocess=NULL))

dev.off()
