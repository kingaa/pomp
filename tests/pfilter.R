options(digits=3)

library(pomp)

ou2() -> ou2

set.seed(9994847L)

pf <- pfilter(ou2,Np=1000)
logLik(pf)
stopifnot(logLik(pf)==sum(cond_logLik(pf)))
range(eff_sample_size(pf))
range(cond_logLik(pf))
names(as(pf,"data.frame"))
identical(as(pf,"data.frame"),as.data.frame(pf))

try(pfilter())
try(pfilter("bob"))

try(pfilter(ou2))
try(pfilter(ou2,Np=NULL))
try(pfilter(ou2,Np=-10))
try(pfilter(ou2,Np=c(10,20,30)))
pfilter(ou2,Np=ceiling(runif(101,min=10,max=100)))

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
stopifnot(all.equal(coef(pf),coef(ou2)))
theta <- coef(ou2)
theta["alpha_2"] <- 0.1
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
theta["tau"] <- 0
try(pfilter(pf,params=theta))
theta <- coef(ou2)
theta["sigma_2"] <- Inf
try(pfilter(pf,params=theta))
theta <- coef(ou2)
theta["alpha_1"] <- 1e60
try(pfilter(pf,params=theta,pred.var=TRUE))

try(pfilter(pf,rprocess=onestep(
  function(x, t, params, delta.t, ...)stop("yikes!"))))
try(pfilter(pf,dmeasure=Csnippet("error(\"ouch!\");")))
pfilter(pf,dmeasure=function(log,...) -Inf)
pfilter(pf,dmeasure=function(log,...) -Inf,filter.mean=TRUE)

pf1 <- pfilter(pf,save.states=TRUE,filter.traj=TRUE)
pf2 <- pfilter(pf,pred.mean=TRUE,pred.var=TRUE,filter.mean=TRUE,save.states="unweighted")
pf3 <- pfilter(pf,t0=1,filter.traj=TRUE)
pf4 <- pfilter(pf,dmeasure=Csnippet("lik = (give_log) ? R_NegInf : 0;"),
  filter.traj=TRUE)
pf1 %>% saved_states() %>% melt() %>% names()
pf1 %>% saved_states() %>% melt() %>% dim()
pf1 %>% saved_states(format="data") %>% names()
pf1 %>% saved_states(format="data") %>% dim()
c(A=pf1,B=pf2) %>% saved_states() %>% melt() %>% names()
c(A=pf1,B=pf2) %>% saved_states() %>% melt() %>% sapply(class)
c(A=pf1,B=pf2) %>% as.data.frame() %>% sapply(class)
c(A=pf1,B=pf2) %>% melt(id=c("time",".id")) %>% sapply(class)
names(as(pf2,"data.frame"))
dim(filter_traj(pf3))
dimnames(filter_traj(pf3))
try(filter_traj(c(pf1,pf3)))
dim(filter_traj(c(pf1,pf4)))
dim(as.data.frame(c(pf1,pf4)))
names(dimnames(filter_traj(c(pf1,pf4))))
names(melt(as(c(pf1,pf4),"data.frame")))
pf2 %>% melt() %>% names()
pf2 %>% melt(id="time") %>% names()
try(saved_states())
try(saved_states(NULL))
try(saved_states("bob"))

try(ou2 %>% as.data.frame() %>% pfilter(Np=1000))

ou2 %>%
  as.data.frame() %>%
  subset(select=c(time,y1,y2)) %>%
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

try(ou2 %>% pfilter(Np=1000) %>% forecast())
try(ou2 %>% pfilter(Np=1000,emeasure=NULL) %>% forecast())
ou2 %>% pfilter(Np=1000,pred.mean=TRUE) %>% forecast() -> y
ou2 %>% pfilter(Np=1000,pred.mean=TRUE) %>% forecast(vars="y1") -> y1
stopifnot(
  dim(y)==c(2,100),
  dim(y1)==c(1,100)
)
