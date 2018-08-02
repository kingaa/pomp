png(filename="probe-%02d.png",res=100)

options(digits=3)
library(pomp)
library(magrittr)

pompExample(gompertz)

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
stopifnot(all.equal(logLik(pb),17.17,tolerance=0.005))

try(probe(po,nsim=10))
try(probe(po,probes=plist[1:3]))
try(probe(po,probes=plist[1:3],nsim=-100))
try(probe(po,probes=plist[1:3],nsim=c(10,20)))
try(probe(po,probes=plist[1:3],nsim=NA))

try(probe(33L))
invisible(probe(pb))
plot(probe(pb,probes=plist[[1]]))
try(probe(pb,probes="okay"))
try(probe(pb,probes=function(x,y)x+y))
try(probe(pb,probes=function(x)stop("hold on now!")))
try({
  count <- 0
  delayed.error <- function (y) {
    count <<- count+1
    if (count>5) stop("whoa nelly!")
    y[1]
  }
  probe(pb,probes=delayed.error)
})
try(probe(pb,nsim=10))

invisible(probe(pb,params=as.list(coef(pb))))
try(probe(pb,params=NULL))
try(probe(pb,params="I think so"))
try({pb1 <- pb; coef(pb1) <- NULL; probe(pb1)})

sapply(probevals(pb),colnames)
sapply(probevals(pb),dim)

dev.off()
