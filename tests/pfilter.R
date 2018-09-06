options(digits=3)

library(pomp2)

pompExample(ou2)

set.seed(9994847L)

pf <- pfilter(ou2,Np=1000)
logLik(pf)
stopifnot(logLik(pf)==sum(cond.logLik(pf)))
range(eff.sample.size(pf))
range(cond.logLik(pf))
names(as(pf,"data.frame"))
identical(as(pf,"data.frame"),as.data.frame(pf))

try(pfilter())
try(pfilter("bob"))

try(pfilter(ou2))
try(pfilter(ou2,Np=NULL))
try(pfilter(ou2,Np=-10))
try(pfilter(ou2,Np=c(10,20,30)))
pfilter(ou2,Np=ceiling(runif(101,min=10,max=100)))
try(pfilter(ou2,Np=100,tol=NA))
try(pfilter(ou2,Np=100,tol=-1))
try(pfilter(ou2,Np=100,tol=NULL))
try(pfilter(ou2,Np=100,tol=c(1,10,20)))

po <- ou2
coef(po) <- NULL
try(pfilter(po,Np=100))
try(pfilter(po,Np=100,params=c()))
try(pfilter(po,Np=100,params=NULL))
try(pfilter(po,Np=100,params=c(1,2,3)))
try(pfilter(po,Np=100,params=c(a=1,b=2,c=3)))
try(pfilter(po,Np=100,params=list()))
pf <- pfilter(po,Np=100,params=as.list(coef(ou2)))

pf <- pfilter(pf)
try(pfilter(pf,Np=-1000))
try(pfilter(pf,tol=-1))
stopifnot(all.equal(coef(pf),coef(ou2)))
theta <- coef(ou2)
theta["alpha.2"] <- 0.1
pf1 <- pfilter(pf,params=theta,Np=100)
stopifnot(identical(coef(pf1),theta))

try(pfilter(ou2,params=parmat(coef(ou2),100)))
try(pfilter(ou2,Np=100,params=parmat(coef(ou2),100)))

pf2 <- pfilter(ou2,Np=function(k)c(100,150)[(k-1)%%2+1])
table(pf2@Np)
try(pfilter(pf2,Np=function(k)c(100,-150)[(k-1)%%2+1]))
try(pfilter(pf2,Np=function(k)c(100,-150)))
try(pfilter(pf2,Np="many"))

theta <- coef(ou2)
theta["sigma.2"] <- Inf
try(pfilter(pf,params=theta))
theta <- coef(ou2)
theta["alpha.1"] <- 1e60
try(pfilter(pf,params=theta,pred.var=TRUE))

try(pfilter(pf,rprocess=onestep(
  function(x, t, params, delta.t, ...)stop("yikes!"))))
try(pfilter(pf,dmeasure=Csnippet("error(\"ouch!\");")))
try(pfilter(pf,dmeasure=function(log,...) -Inf))

set.seed(388966382L)
capture.output(try(pfilter(pf,Np=2,max.fail=20,verbose=TRUE,filter.mean=TRUE)),
  type="message") -> out
stopifnot(sum(grepl("filtering failure at",out))==21)
stopifnot(grepl("too many filtering failures",tail(out,1)))

pf1 <- pfilter(pf,save.states=TRUE,filter.traj=TRUE)
pf2 <- pfilter(pf,pred.mean=TRUE,pred.var=TRUE,filter.mean=TRUE)
pf3 <- pfilter(pf,t0=1,filter.traj=TRUE)
pf4 <- pfilter(pf,dmeasure=Csnippet("lik = (give_log) ? R_NegInf : 0;"),
  filter.traj=TRUE)
names(as(pf2,"data.frame"))
dim(filter.traj(pf3))
dimnames(filter.traj(pf3))
try(filter.traj(c(pf1,pf3)))
dim(filter.traj(c(pf1,pf4)))
names(dimnames(filter.traj(c(pf1,pf4))))

library(magrittr)

try(ou2 %>% as.data.frame() %>% pfilter(Np=1000))

ou2 %>%
  as.data.frame() %>%
  pfilter(
    times="time",t0=0,Np=500,
    params=list(x1_0=-3,x2_0=4),
    rprocess=onestep(
    step.fun=function(x1,x2,delta.t,...) {
      setNames(rnorm(n=2,mean=c(x1,x2),sd=5*delta.t),c("x1","x2"))
    }
  ),
    dmeasure=function(x1,x2,y1,y2,...,log) {
      ll <- sum(dnorm(x=c(y1,y2),mean=c(x1,x2),sd=5,log=TRUE))
      if (log) ll else exp(ll)
    }
  )

