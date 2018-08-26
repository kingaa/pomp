options(digits=3)

library(pomp)
library(magrittr)

pompExample(ou2,envir=NULL) -> ou2
ou2[[1]] -> po

coef(po,"alpha.sd") <- 5

set.seed(1835425749L)

po %>%
  pomp(
    dprior=function(alpha.1,alpha.2,alpha.3,alpha.4,alpha.sd,...,log) {
      ll <- sum(
        dnorm(
          x=c(alpha.1,alpha.2,alpha.3,alpha.4),
          mean=c(0.8,-0.5,0.3,0.9),
          sd=alpha.sd,
          log=TRUE
        )
      )
      if (log) ll else exp(ll)
    },
    rprior=function(params,...) {
      params[c("alpha.1","alpha.2","alpha.3","alpha.4")] <- rnorm(
        n=4,
        mean=c(0.8,-0.5,0.3,0.9),
        sd=params["alpha.sd"]
      )
      params
    }
  ) -> po

stopifnot(
  po %>% dprior(params=coef(po),log=TRUE) == 4*dnorm(x=0,mean=0,sd=5,log=TRUE),
  all.equal(po %>% dprior(params=coef(po)),dnorm(x=0,mean=0,sd=5)^4)
)

replicate(5,rprior(po,params=coef(po))) %>% parmat() -> theta
stopifnot(round(dprior(po,params=theta,log=TRUE),3) ==
  c(-12.237, -10.848, -15.806, -10.847, -11.526))
