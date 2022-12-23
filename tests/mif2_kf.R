options(digits=3)

library(pomp)
suppressPackageStartupMessages({
  library(dplyr)
})

set.seed(376098756)

gompertz() |> window(end=10) |> simulate(seed=1176423047) -> po

po |>
  mif2(Nmif=100,Np=1000,cooling.fraction.50=0.4,cooling.type="geometric",
    rw.sd=rw_sd(sigma=0.02,r=0.02,X_0=ivp(0.05),tau=0.02)) |>
  continue(Nmif=100) -> mf
replicate(n=10,mf |> pfilter(Np=3000)) -> pfs
pfs |> sapply(logLik) |> logmeanexp(se=TRUE) -> pf.ll.mle

replicate(n=10,po |> pfilter(Np=3000)) |>
  sapply(logLik) |>
  logmeanexp(se=TRUE) -> pf.ll.truth

po |>
  as.data.frame() |>
  mutate(logY=log(Y)) |>
  select(time,logY) |>
  pomp(t0=0,times="time") -> logpo
  
po |> coef() |> as.list() -> theta
po |> rinit(params=coef(po)) |> as.numeric() -> x0
with(theta,c(x=log(x0/K))) -> X0
po |> obs() |> log() -> y
with(theta,matrix(c(exp(-r)),1,1)) -> A
with(theta,matrix(c(sigma*sigma),1,1)) -> Q
with(theta,matrix(1,1,1)) -> C
with(theta,tau*tau) -> R
kalmanFilter(logpo,X0=X0,A=A,Q=Q,C=C,R=R) -> kf.truth

mf |> coef() |> as.list() -> theta
mf |> rinit(params=coef(mf)) |> as.numeric() -> x0
mf |> obs() |> log() -> y
with(theta,c(x=log(x0/K))) -> X0
with(theta,matrix(c(exp(-r)),1,1)) -> A
with(theta,matrix(c(sigma*sigma),1,1)) -> Q
with(theta,matrix(1,1,1)) -> C
with(theta,tau*tau) -> R
kalmanFilter(logpo,X0=X0,A=A,Q=Q,C=C,R=R) -> kf.mle

cat("likelihood at truth:",kf.truth$logLik-sum(y),"\n")
cat("pfilter likelihood at truth:",pf.ll.truth,"\n")
cat("likelihood at IF2 mle:",kf.mle$logLik-sum(y),"\n")
cat("pfilter likelihood at IF2 mle:",pf.ll.mle,"\n")
