### OU2 test of abc for pomp

png(filename="abc-%02d.png",res=100)

library(pomp)

ou2() -> ou2

set.seed(2079015564L)

plist <- list(
  y1.mean=probe.mean(var="y1"),
  y2.mean=probe.mean(var="y2"),
  probe.acf(var="y1",lags=c(0,5)),
  probe.acf(var="y2",lags=c(0,5)),
  probe.ccf(vars=c("y1","y2"),lags=0)
)

ou2 %>% probe(probes=plist,nsim=100) -> pb

sqrt(diag(covmat(pb))) -> scale.dat

ou2 %>%
  abc(Nabc=100,probes=plist,scale=scale.dat,epsilon=1.7,
    proposal=mvn.diag.rw(rw.sd=c(alpha_1=0.01,alpha_2=0.01)),
    dprior=function(alpha_1,alpha_2,...,log) {
      ll <- sum(dnorm(x=c(alpha_1,alpha_2),mean=c(0.6,0),sd=4,log=TRUE))
      if (log) ll else exp(ll)
    }) -> abc1
plot(abc1)
plot(abc1,scatter=TRUE)

crossprod(
  array(data=c(0.1,0.02,0,0.1),dim=c(2,2),
    dimnames=list(c("alpha_1","alpha_2"),c("alpha_1","alpha_2")))
) -> sig

pb %>% abc(Nabc=100,scale=scale.dat,epsilon=2,proposal=mvn.rw(sig)) -> abc2
abc2 %>% abc(Nabc=100) -> abc3
abc1 %>% abc(Nabc=80) %>% continue(Nabc=20) -> abc4

plot(c(abc1,abc2,abc3,abc4),y="bob")
plot(c(abc1,abc2,abc3,abc4),scatter=TRUE)

c(a=c(abc1,abc2),b=abc3) -> abclist
stopifnot(identical(abclist,c(a1=abc1,c(a2=abc2,b=abc3))))
stopifnot(all(dim(traces(abc1))==c(101,10)))
stopifnot(all(dim(traces(abc1,"alpha_1"))==c(101,1)))
invisible(conv.rec(abc2))
dim(as.data.frame(abclist))
  
c(abc1,abc2) %>% traces() -> traces
traces %>% length()
traces %>% class()
traces %>% sapply(dim)
try(abclist %>% plot(pars="alpha_3",scatter=TRUE))

abc1 %>%
  abc(Nabc=500,dprior=Csnippet("
    lik = dnorm(alpha_1,0.8,1,1)+dnorm(alpha_2,0.2,1,1);
    lik = (give_log) ? lik : exp(lik);"
  ),paramnames=c("alpha_1","alpha_2")) -> abc5

abc1 %>% abc(Nabc=50,params=as.list(coef(ou2))) %>% plot()

abc4 %>% abc(proposal=function(theta,...)theta) %>% plot()

try(abc())
try(abc(3))

s5 <- simulate(abc5)
stopifnot(
  is(abc5,"abcd_pomp"),
  is(simulate(s5),"pomp"),
  !is(simulate(s5),"abcd_pomp")
)

try(abc(abc1,Nabc=-5))
stopifnot(all(dim(traces(abc(abc1,Nabc=0))==c(1,10))))

try(abc(ou2,Nabc=50,scale=scale.dat[1:2],epsilon=1.7,
  proposal=mvn.diag.rw(rw.sd=c(alpha_1=0.01,alpha_2=0.01))))
try(abc(ou2,Nabc=50,probes=plist,scale=scale.dat[1:2],epsilon=1.7,
  proposal=mvn.diag.rw(rw.sd=c(alpha_1=0.01,alpha_2=0.01))))
po <- ou2
coef(po) <- NULL
try(abc(po,Nabc=100,probes=plist,scale=scale.dat,epsilon=1.7,
  proposal=mvn.diag.rw(rw.sd=c(alpha_1=0.01,alpha_2=0.01))))
try(abc(ou2,params=numeric(0),Nabc=100,probes=plist,scale=scale.dat,
  epsilon=1.7,proposal="mvn.diag.rw(rw.sd=c(alpha_1=0.01,alpha_2=0.01))"))
try(abc(ou2,params=NULL,Nabc=100,probes=plist,scale=scale.dat,epsilon=1.7,
  proposal=mvn.diag.rw(rw.sd=c(alpha_1=0.01,alpha_2=0.01))))
try(abc(ou2,Nabc=100,probes=plist,scale=scale.dat,epsilon=1.7,
  proposal="mvn.diag.rw(rw.sd=c(alpha_1=0.01,alpha_2=0.01))"))
try(abc(ou2,Nabc=100,probes=plist,scale=scale.dat,epsilon=1.7,
  proposal=function(...)stop("yikes!")))
try(abc(ou2,Nabc=100,probes=plist,scale=scale.dat,epsilon=1.7,
  proposal=function(...)3))
try(abc(ou2,Nabc=100,probes=plist,scale=scale.dat,epsilon=1.7))
try(abc(ou2,Nabc=100,proposal=function(theta,...)theta,probes="mary",
  scale=scale.dat,epsilon=1.7))
try(abc(ou2,Nabc=100,proposal="bob",probes="mary",epsilon=1.7))
try(abc(ou2,Nabc=100,probes=plist,epsilon=1.7,
  proposal=mvn.diag.rw(rw.sd=c(alpha_1=0.01,alpha_2=0.01))))
try(abc(ou2,Nabc=100,probes=plist,scale=scale.dat,
  proposal=mvn.diag.rw(rw.sd=c(alpha_1=0.01,alpha_2=0.01))))
try(abc(ou2,Nabc=100,probes=plist,scale=scale.dat,
  proposal=mvn.diag.rw(rw.sd=c(alpha_1=0.01,alpha_2=0.01)),
  epsilon=1,rprocess=NULL))

try(abc(abc1,Nabc=100,epsilon=NULL,scale=scale.dat))
try(abc(ou2,params=c(1,2,3),Nabc=100,probes=plist,scale=scale.dat,epsilon=1.7,
  proposal=mvn.diag.rw(rw.sd=c(alpha_1=0.01,alpha_2=0.01))))
try(abc(ou2,Nabc=100,probes="plist[[1]]",scale=scale.dat[1],epsilon=1.7,
  proposal=mvn.diag.rw(rw.sd=c(alpha_1=0.01,alpha_2=0.01))))
try(abc(ou2,Nabc=100,probes=function(x,y)x+y,scale=scale.dat[1],epsilon=1.7,
  proposal=mvn.diag.rw(rw.sd=c(alpha_1=0.01,alpha_2=0.01))))
try(abc(abc1,dprior=function(log,...)stop("ouch!")))
try(abc(abc1,dprior=function(log,...)Inf))
try(abc(abc1,probes=function(x)stop("piff!")))
count <- 0
delayed.failure <- function (x) {count <<- count+1; if (count>2) stop("paff!") else 1}
try(abc(abc1,scale=1,probes=delayed.failure))
try(abc(abc1,proposal=function(...)stop("urp!")))
count <- 0
delayed.failure <- function (log,...) {count <<- count+1; if (count>5) stop("no sir!") else 1}
try(abc(abc1,dprior=delayed.failure))
count <- 0
delayed.failure <- function (theta,...) {count <<- count+1; if (count>5) stop("'fraid not!") else theta}
try(abc(abc1,proposal=delayed.failure))

coef(c(abc1,ou2)) -> theta
stopifnot(dim(theta)==c(10,2),
  rownames(theta)==c("alpha_1","alpha_2","alpha_3","alpha_4",
    "sigma_1","sigma_2","sigma_3","tau","x1_0","x2_0"),
  apply(theta[c("alpha_3","alpha_4",
    "sigma_1","sigma_2","sigma_3","tau","x1_0","x2_0"),],1,diff)==0,
  apply(theta[c("alpha_1","alpha_2"),],1,diff) != 0)
try(c(abc1,NULL))
try(c(c(abc1,abc2),ou2))
c(abc1)
alist <- c(c(abc1,abc2))
class(alist[2])
try(alist[3])
alist <- c(a=abc1,b=abc2)
alist["b"]
alist["c"]
alist[["b"]]
alist[["c"]]
c(one=abc1,two=abc2,three=abc3)
print(c(one=abc1,two=abc2,three=abc3))

capture.output(abc(ou2,Nabc=100,probes=plist,scale=scale.dat,epsilon=1.7,
  proposal=mvn.diag.rw(rw.sd=c(alpha_1=0.01,alpha_2=0.01)),
  verbose=TRUE) -> abc1) -> out
stopifnot(
  length(out)==40,
  sum(grepl("acceptance",out))==20,
  sum(grepl("ABC iteration",out))==20
)

gompertz() -> gompertz
set.seed(2079015564L)

gompertz %>%
  as.data.frame() %>%
  abc(Nabc=20,times="time",t0=0,
    scale=1,epsilon=10,
    probes=list(probe.mean("Y"),probe.median("Y")),
    partrans=parameter_trans(log=c("r","K")),
    paramnames=c("r","K"),
    proposal=mvn.diag.rw(rw.sd=c(r=0.01,K=0.01)),
    params=coef(gompertz),
    rinit=function(...)c(X=1),
    rprocess=discrete_time(function (X, r, K, ...) c(X=r*X*exp(-X/K))),
    rmeasure=function (Y, X, ...) c(Y = rnorm(n=1,mean=X,sd=2))
    ) %>% plot()

dev.off()
