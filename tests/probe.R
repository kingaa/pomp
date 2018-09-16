options(digits=3)
png(filename="probe-%02d.png",res=100)

library(pomp2)
library(magrittr)

gompertz() -> gompertz

set.seed(234501470L)

po <- gompertz

plist <- list(
  mean=probe.mean("Y",trim=0.1,transform=sqrt),
  sd=probe.sd("Y",transform=sqrt),
  probe.marginal("Y",ref=obs(po)),
  probe.acf("Y",lags=c(1,3,5),type="correlation",transform=sqrt),
  probe.quantile("Y",prob=c(0.25,0.75))
)

probe(po,probes=plist,nsim=500,seed=595969) -> pb
plot(pb,y=NULL)
pb %>% values() %>% head(3) %>% knitr::kable()
pb %>% as.data.frame() %>% head(3) %>% knitr::kable()
summary(pb)

try(probe())
try(probe("po"))
try(probe(NULL))
try(probe(po,nsim=10))
try(probe(po,probes=plist[1:3]))
try(probe(po,probes=plist[1:3],nsim=-100))
try(probe(po,probes=plist[1:3],nsim=c(10,20)))
try(probe(po,probes=plist[1:3],nsim=NA))
try(probe(po,nsim=100,probes=function(x)rep(1,times=ceiling(runif(1,max=10)))))

try(probe(33L))
probe(pb)
plot(probe(pb,probes=plist[[1]]))
try(probe(pb,probes="okay"))
try(probe(pb,probes=function(x,y)x+y))
try(probe(pb,probes=function(x)stop("hold on now!")))
try(probe(pb,probes=function(x)"bob"))
try({
  count <- 0
  delayed.error <- function (y) {
    count <<- count+1
    if (count>5) stop("whoa nelly!")
    y[1]
  }
  probe(pb,probes=delayed.error)
})
try({
  count <- 0
  delayed.badval <- function (y) {
    count <<- count+1
    if (count>10) "bob" else 10
  }
  probe(pb,probes=delayed.badval)
})
try({
  count <- 0
  delayed.badval <- function (y) {
    count <<- count+1
    if (count>10) rep(3.5,count) else 3.5
  }
  probe(pb,probes=delayed.badval)
})

try(probe(pb,nsim=10))

probe(pb,params=as.list(coef(pb)))
try(probe(pb,params=NULL))
try(probe(pb,params="I think so"))
try({pb1 <- pb; coef(pb1) <- NULL; probe(pb1)})

po %>% probe(nsim=100,probes=function(x)1) %>% logLik()

try(data.frame(t=1:10,a=1:10) %>% probe())

data.frame(t=1:10,a=1:10) %>%
  probe(
    times="t",t0=0,
    rprocess=euler(
      function(t,x,delta.t,...){
        c(x=rlnorm(n=1,meanlog=log(x),sd=sqrt(delta.t)))
      },delta.t=0.1),
    rmeasure=function(t,x,...){
      c(a=rpois(n=1,lambda=x))
    },
    nsim=1000,
    params=c(x_0=1),
    probes=list(
      f=probe.mean("a",transform=sqrt),
      g=probe.median("a"),
      h=function(y)range(y)
    )
  )

ou2() -> ou2
ou2 %>% probe(nsim=100,probes=probe.ccf(c("y1","y2"),lags=c(-10,0,1))) %>% plot()

dev.off()
