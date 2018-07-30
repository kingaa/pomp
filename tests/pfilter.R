options(digits=4)

library(pomp)

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

pp <- parmat(coef(ou2),100)
pf <- pfilter(ou2,params=pp)
pf <- pfilter(pf1,params=pp)
try(pfilter(ou2,params=pp,Np=1000))
rownames(pp) <- NULL
try(pfilter(ou2,params=pp))

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

try(pfilter(pf,rprocess=onestep.sim(
  function(x, t, params, delta.t, ...)stop("yikes!"))))
try(pfilter(pf,dmeasure=Csnippet("error(\"ouch!\");")))
try(pfilter(pf,dmeasure=function(y,x,t,params,log,...) -Inf))

set.seed(388966382L)
capture.output(try(pfilter(pf,Np=2,max.fail=20,verbose=TRUE)),
  type="message") -> out
stopifnot(sum(grepl("filtering failure at",out))==21)
stopifnot(grepl("too many filtering failures",tail(out,1)))

pf1 <- pfilter(pf,save.states=TRUE,filter.traj=TRUE,save.params=TRUE)
pf2 <- pfilter(pf,pred.mean=TRUE,pred.var=TRUE,filter.mean=TRUE)
pf3 <- pfilter(pf,t0=1,filter.traj=TRUE)
pf4 <- pfilter(pf,dmeasure=Csnippet("lik = (give_log) ? R_NegInf : 0;"),
               filter.traj=TRUE)
names(as(pf2,"data.frame"))
dim(filter.traj(pf3))
dimnames(filter.traj(pf3))
