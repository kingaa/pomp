options(digits=3)

library(pomp)
suppressPackageStartupMessages({
  library(dplyr)
})

ou2() -> po

stopifnot(
  po |> rprior() == coef(po)
)

coef(po,"alpha_sd") <- 5

set.seed(1835425749L)

po |>
  pomp(
    dprior=function(alpha_1,alpha_2,alpha_3,alpha_4,alpha_sd,...,log) {
      ll <- sum(
        dnorm(
          x=c(alpha_1,alpha_2,alpha_3,alpha_4),
          mean=c(0.8,-0.5,0.3,0.9),
          sd=alpha_sd,
          log=TRUE
        )
      )
      if (log) ll else exp(ll)
    },
    rprior=function(alpha_1,alpha_2,alpha_3,alpha_4,alpha_sd,...) {
      c(
        alpha_1=rnorm(n=1,mean=0.8,sd=alpha_sd),
        alpha_2=rnorm(n=1,mean=-0.5,sd=alpha_sd),
        alpha_3=rnorm(n=1,mean=0.3,sd=alpha_sd),
        alpha_4=rnorm(n=1,mean=0.9,sd=alpha_sd)
      )
    }
  ) -> po

stopifnot(
  po |> dprior(params=coef(po),log=TRUE) == 4*dnorm(x=0,mean=0,sd=5,log=TRUE),
  all.equal(po |> dprior(params=coef(po)),dnorm(x=0,mean=0,sd=5)^4)
)

replicate(5,rprior(po,params=coef(po))) |> parmat() -> theta
stopifnot(round(dprior(po,params=theta,log=TRUE),3) ==
            c(-12.237, -10.848, -15.806, -10.847, -11.526))

try(dprior("ou2",params=theta))
try(dprior(params=theta))

try(rprior("ou2",params=theta))
try(rprior(params=theta))

try(po |> pomp(rprior=function(...)c(1,3,3)) |> rprior(params=theta))
po |> rprior(params=parmat(theta,3)) -> p
stopifnot(
  dim(p)==c(11,15),
  names(dimnames(p))==c("name",".id"),
  rownames(p)==names(coef(po))
)

freeze(
  mean(dprior(po,params=replicate(n=10000,rprior(po)),log=TRUE)),
  seed=1865498175
) -> L
stopifnot(round(L,1)==-12.1)
