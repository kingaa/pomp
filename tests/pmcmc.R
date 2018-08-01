options(digits=3)
png(filename="pmcmc-%02d.png",res=100)

set.seed(857075216L)

library(pomp)

pompExample(gompertz,envir=NULL) -> po

po <- window(po[[1]],end=10)

prop1 <- mvn.diag.rw(c(r=0.01,sigma=0.01))

mcmc1 <- pmcmc(po,Nmcmc=100,Np=100,dprior=Csnippet("
    lik = dunif(r,0,1,1)+dnorm(sigma,0,1,1);
    lik = (give_log) ? lik : exp(lik);"
),
  paramnames=c("r","sigma"),
  proposal=prop1)

try(pmcmc(po,Nmcmc=10,Np="bob"))
try(pmcmc(po,Nmcmc=10))
try(pmcmc(po,Nmcmc=10,Np=NULL))
try(pmcmc(po,Nmcmc=10,Np=function(k)10*k))
try(pmcmc(po,Nmcmc=10,Np=function(k)NULL))
try(pmcmc(po,Nmcmc=10,Np=rep(100,20)))
try(pmcmc(po,Nmcmc=10,Np=-50))
try(pmcmc(po,Nmcmc=10,Np=NA))
try(pmcmc(po,Nmcmc=10,Np=function(k)10*(k+1)))
try(pmcmc(po,Nmcmc=100,Np=100))

pf <- pfilter(po,Np=100)
mcmc2 <- pmcmc(pf,Nmcmc=100,proposal=prop1)
mcmc3 <- pmcmc(mcmc1,Nmcmc=50)
mcmc3 <- continue(mcmc3,Nmcmc=50)

plot(c(mcmc1,mcmc2,mcmc3),pars=c("r","sigma"),density=FALSE)
plot(c(mcmc1,c(mcmc2,mcmc3)),pars=c("r","sigma"),trace=FALSE)
window(traces(c(c(mcmc1,c(mcmc2,mcmc3)))),thin=50)
plot(traces(c(c(mcmc1,mcmc2),mcmc3),c("r","sigma")))
try(plot(traces(c(c(mcmc1,mcmc2),mcmc3),c("r","bob"))))
print(mcmc1)
c(mcmc1,mcmc2) -> mcl
mcl[1]
mcl[3]

try(pmcmc())
try(pmcmc("bob"))

invisible(pmcmc(mcmc1,start=as.list(coef(mcmc3))))
try(pmcmc(mcmc1,start=NULL))
try(pmcmc(mcmc1,start=-7))
try(pmcmc(mcmc1,start="yes"))
try(pmcmc(mcmc1,start=list()))
try({tmp <- mcmc1; coef(tmp) <- NULL; pmcmc(tmp)})
try(pmcmc(mcmc1,proposal="random"))
try(pmcmc(mcmc1,proposal=NULL))
try(pmcmc(mcmc1,proposal=function(...)c(3,2)))
try(pmcmc(mcmc1,proposal=function(...)c(a=3,b=2,X.0=1)))
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

try(pmcmc(mcmc1,dprior=function(params,log,...)stop("not again!")))
try({
  count <- 0
  delayed.failure <- function (params, log, ...) {
    count <<- count+1
    if (count>5) stop("uh huh") else 1
  }
  pmcmc(mcmc1,dprior=delayed.failure)})

capture.output(invisible(pmcmc(mcmc1,Nmcmc=10,verbose=TRUE))) -> out
stopifnot(sum(grepl("acceptance ratio",out))==10)
stopifnot(sum(grepl("PMCMC iteration",out))==11)

dev.off()
