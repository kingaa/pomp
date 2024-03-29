options(digits=3)
png(filename="pmcmc-%02d.png",res=100)

library(pomp)

set.seed(857075216L)

gompertz() |> window(end=10) -> po

prop1 <- mvn_diag_rw(c(r=0.01,sigma=0.01))

mcmc1 <- pmcmc(
  po,Nmcmc=100,Np=100,dprior=Csnippet("
    lik = dunif(r,0,1,1)+dnorm(sigma,0,1,1);
    lik = (give_log) ? lik : exp(lik);"),
  paramnames=c("r","sigma"),
  proposal=prop1
)

try(pmcmc())
try(pmcmc("bob"))
try(pmcmc(po))
try(pmcmc(po,proposal="yes"))
try(pmcmc(po,proposal=NULL))
try(pmcmc(po,proposal=3))
try(pmcmc(po,proposal="prop1"))
try(pmcmc(po,proposal=prop1))
try(pmcmc(po,proposal=prop1,Np="bob"))
try(pmcmc(po,proposal=prop1,Np=NULL))
try(pmcmc(po,proposal=prop1,Np=NA))
try(pmcmc(po,proposal=prop1,Np=c(1,2,3)))
try(pmcmc(po,proposal=prop1,Np=-5))
pmcmc(po,proposal=prop1,Np=1)
pmcmc(po,proposal=prop1,Np=rep(1,11))
try(pmcmc(po,proposal=prop1,Np=function(k)10*k))
pmcmc(po,proposal=prop1,Np=function(k)10*k+10)

try(pmcmc(mcmc1,params=c(r=3,K=1,sigma=0.1,tau=0.1,X_0=1),Np=100))

pf <- pfilter(po,Np=100)
mcmc2 <- pmcmc(pf,Nmcmc=100,proposal=prop1)
mcmc3 <- pmcmc(mcmc1,Nmcmc=50)
mcmc3 <- continue(mcmc3,Nmcmc=50)

plot(c(mcmc1,mcmc2,mcmc3),pars=c("r","sigma"),density=FALSE)
plot(c(mcmc1,c(mcmc2,mcmc3)),pars=c("r","sigma"),trace=FALSE)
invisible(window(traces(c(c(mcmc1,c(mcmc2,mcmc3)))),thin=10))
plot(traces(c(c(mcmc1,mcmc2),mcmc3),c("r","sigma")))
try(traces(c(c(mcmc1,mcmc2),mcmc3),c("r","bob")))
dim(as.data.frame(c(mcmc1,mcmc2)))

filter_traj(c(mcmc1,mcmc2,mcmc3)) -> ft
stopifnot(
  dim(ft)==c(1,100,11,3),
  names(dimnames(ft))==c("name","rep","time","chain")
)

print(mcmc1)
c(mcmc1,mcmc2) -> mcl
mcl[1]
mcl[3]

stopifnot(dim(filter_traj(mcmc1))==c(1,100,11),
  dim(filter_traj(c(mcmc1,mcmc2,mcmc3)))==c(1,100,11,3))
logLik(mcmc1)
logLik(c(mcmc1,mcmc2,mcmc3))

pmcmc(mcmc1,params=as.list(coef(mcmc3)))
try(pmcmc(mcmc1,params=NULL))
try(pmcmc(mcmc1,params=-7))
try(pmcmc(mcmc1,params="yes"))
try(pmcmc(mcmc1,params=list()))
try({tmp <- mcmc1; coef(tmp) <- NULL; pmcmc(tmp)})
try(pmcmc(mcmc1,proposal=function(...)c(3,2)))
try(pmcmc(mcmc1,Nmcmc=3,proposal=function(...)c(a=3,b=2,X.0=1)))
try(pmcmc(mcmc1,proposal=function(...)stop("oh no!")))
try({
  count <- 0
  delayed.failure <- function (theta, ...) {
    count <<- count+1
    if (count>5) stop("no sir!") else theta
  }
  pmcmc(mcmc1,proposal=delayed.failure)})

try(pmcmc(mcmc1,Nmcmc=-20))
try(pmcmc(mcmc1,Nmcmc=NA))
try(pmcmc(mcmc1,Nmcmc=c(5,10,15)))

try(pmcmc(mcmc1,dprior=function(log,...)stop("not again!")))
try({
  count <- 0
  delayed.failure <- function (log, ...) {
    count <<- count+1
    if (count>5) stop("uh huh") else 1
  }
  pmcmc(mcmc1,dprior=delayed.failure)})

try(pmcmc(mcmc1,dmeasure=function(log,...)if (log) -Inf else 0))

capture.output(invisible(pmcmc(mcmc1,Nmcmc=10,verbose=TRUE))) -> out
stopifnot(sum(grepl("acceptance ratio",out))==10)
stopifnot(sum(grepl("PMCMC iteration",out))==11)

gompertz() -> gompertz

set.seed(857075216L)

try(gompertz |> as.data.frame() |> pmcmc())
try(gompertz |> as.data.frame() |> pmcmc(times="time",t0=0))
try(gompertz |> as.data.frame() |> pmcmc(times="time",t0=0,
  proposal=mvn_diag_rw(c(a=1,b=2))))

gompertz |>
  as.data.frame() |>
  pmcmc(
    Nmcmc=10,Np=100,
    times="time",t0=0,
    rprocess=discrete_time(
      function (x, r, K, ...) {
        c(x=x*exp(r*(1-x/K)))
      }
    ),
    dmeasure=function (Y, x, ..., log) {
      dlnorm(Y,meanlog=log(0.01*x),sdlog=2,log=log)
    },
    dprior=function(r,K,...,log) {
      ll <- sum(dlnorm(x=c(r,K),meanlog=log(0.1,150),sdlog=3,log=TRUE))
      if (log) ll else exp(ll)
    },
    proposal=mvn_diag_rw(c(r=0.01,K=10)),
    params=c(r=0.1,K=150,x_0=150)
  ) -> mcmc5
plot(mcmc5,pars=c("r","K"))

dev.off()
