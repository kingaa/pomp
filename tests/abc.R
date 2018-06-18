### OU2 test of abc for pomp

library(pomp)
pompExample(ou2)

pdf(file='abc.pdf')

set.seed(2079015564L)

probes.good <- list(
    y1.mean=probe.mean(var="y1"),
    y2.mean=probe.mean(var="y2"),
    probe.acf(var="y1",lags=c(0,5)),
    probe.acf(var="y2",lags=c(0,5)),
    probe.ccf(vars=c("y1","y2"),lags=0)
)
psim <- probe(ou2,probes=probes.good,nsim=200)
plot(psim)
## why do simulations sometimes seem extreme with respect to these probes?

scale.dat <- apply(psim$simvals,2,sd)

po <- ou2

abc1 <- abc(po,
            Nabc=2000,
            start=coef(ou2),
            probes=probes.good,
            scale=scale.dat,
            epsilon=1.7,
            proposal=mvn.diag.rw(rw.sd=c(alpha.1=0.01,alpha.2=0.01))
            )

plot(abc1,scatter=TRUE)
plot(abc1)
plot(abc1,pars=c("alpha.1","alpha.3"))
abc1@pars <- character(0)
plot(abc1)

try(abc(po,Nabc=2000,probes=probes.good,start=numeric(0),scale=scale.dat,
        epsilon=1.7,proposal=mvn.diag.rw(rw.sd=c(alpha.1=0.01,alpha.2=0.01))))
try(abc(po,Nabc=2000,probes=probes.good,start=unname(coef(po)),scale=scale.dat,
  epsilon=1.7,proposal=mvn.diag.rw(rw.sd=c(alpha.1=0.01,alpha.2=0.01))))
try(abc(po,Nabc=2000,probes=probes.good,start=numeric(0),scale=scale.dat,
  epsilon=1.7,proposal=mvn.diag.rw(rw.sd=c(alpha.1=0.01,alpha.2=0.01))))
try(abc(po,Nabc=2000,probes=probes.good,scale=scale.dat,epsilon=1.7,proposal="bob"))
try(abc(po,Nabc=2000,probes="probes.good",scale=scale.dat,
        epsilon=1.7,proposal=mvn.diag.rw(rw.sd=c(alpha.1=0.01,alpha.2=0.01))))
try(abc(po,Nabc=2000,probes=function(x,y)x,scale=scale.dat,
        epsilon=1.7,proposal=mvn.diag.rw(rw.sd=c(alpha.1=0.01,alpha.2=0.01))))
try(abc(pomp(po,dprior=function(params,log,...)Inf),Nabc=2000,probes=probes.good,
        scale=scale.dat,epsilon=1.7,
        proposal=mvn.diag.rw(rw.sd=c(alpha.1=0.01,alpha.2=0.01))))
try(abc(pomp(po,Nabc=2000,probes=probes.good,
        scale=scale.dat,epsilon=1.7)))
try(abc(po,Nabc=2000,scale=scale.dat,
        epsilon=1.7,proposal=mvn.diag.rw(rw.sd=c(alpha.1=0.01,alpha.2=0.01))))
try(abc(po,Nabc=2000,probes=probes.good,
        epsilon=1.7,proposal=mvn.diag.rw(rw.sd=c(alpha.1=0.01,alpha.2=0.01))))
try(abc(po,Nabc=2000,probes=probes.good,scale=scale.dat,
        proposal=mvn.diag.rw(rw.sd=c(alpha.1=0.01,alpha.2=0.01))))


## check how sticky the chain is:
runs <- rle(as.vector(conv.rec(abc1)[, "alpha.1"]))
hist(runs$lengths)
mean(runs$length)

abc2 <- abc(po,
            Nabc=500,
            probes=probes.good,
            scale=scale.dat,
            epsilon=1,
            proposal=mvn.diag.rw(c(alpha.1=0.01,alpha.2=0.01))
            )
plot(abc2)

abc3 <- abc(po,
            Nabc=500,
            probes=probes.good,
            scale=scale.dat,
            epsilon=2,
            proposal=mvn.diag.rw(c(alpha.1=0.01,alpha.2=0.01))
            )
abc3 <- continue(abc3,Nabc=700)
plot(abc3)

sig <- array(data=c(0.1,0.02,0,0.1),
             dim=c(2,2),
             dimnames=list(c("alpha.1","alpha.2"),
                           c("alpha.1","alpha.2")))
sig <- crossprod(sig)

abc4 <- abc(probe(po,probes=probes.good,nsim=200),
            Nabc=500,
            scale=scale.dat,
            epsilon=2,
            proposal=mvn.rw(sig)
            )
plot(abc4)

abc5 <- abc(abc4,Nabc=250)
plot(abc5)
abc5 <- abc(abc5)

dprior6 <- function (params, log, ...) {
    ll <- sum(
        dnorm(
            x=params[c("alpha.1","alpha.2","alpha.3","alpha.4")],
            mean=c(0.8,-0.5,0.3,0.9),
            sd=5,
            log=TRUE
        )
    )
    if (log) ll else exp(ll)
}

abc6 <- abc(pomp(po,dprior=dprior6),
            Nabc=500,
            probes=probes.good,
            scale=scale.dat,
            epsilon=1,
            proposal=mvn.diag.rw(c(alpha.1=0.01,alpha.2=0.01))
            )
plot(abc6)

try(c(abc6,ou2))
try(c(c(abc6,abc6),ou2))
abc7 <- c(abc2)
abc7 <- c(abc7)
try(abc7 <- c(abc2,abc3))
try(abc7 <- c(abc7,abc3))
plot(abc7 <- c(abc2,abc4))
abc2 <- abc7[1]
abc2 <- abc2[[1]]
plot(abc2,scatter=TRUE,y=NA)
plot(abc7,scatter=TRUE,y=NA)
try(plot(abc7,pars=c("alpha.1"),scatter=TRUE,y=NA))
plot(conv.rec(c(abc2,abc4)))
plot(conv.rec(c(abc7,abc6)))
plot(window(conv.rec(c(abc7,abc6),c("alpha.1","alpha.2")),thin=20,start=100))
invisible(covmat(abc7))

capture.output(
    abc8 <- abc(
        pomp(ou2,dprior=function (params, log, ...) {
            f <- sum(dnorm(params,mean=coef(ou2),sd=1,log=TRUE))
            if (log) f else exp(f)
        }),
        Nabc=500,verbose=TRUE,
        probes=probes.good,
        scale=scale.dat,
        epsilon=5,
        proposal=mvn.rw.adaptive(rw.sd=c(alpha.2=0.01,alpha.3=0.01),
                                 scale.start=500,shape.start=100))
) -> out
stopifnot(length(out)==2201)
stopifnot(sum(grepl("acceptance ratio",out))==100)

abc8 <- continue(abc8,Nabc=2000,proposal=mvn.rw(covmat(abc8)))

plot(abc8,scatter=TRUE)
plot(abc8)

traces <- window(conv.rec(abc8,c("alpha.2","alpha.3")),start=500)
library(coda)
rejectionRate(traces)
autocorr.diag(traces)
traces <- window(traces,thin=50)
geweke.diag(traces)

vmat <- matrix(c(1,0,0,0),ncol=2,nrow=2,
               dimnames=list(c("alpha.1","alpha.2"),
                             c("alpha.1","alpha.2")))
abc9 <- abc(
    abc8,Nabc=1,probes=probes.good,scale=scale.dat,epsilon=5,
    proposal=mvn.rw.adaptive(rw.var=vmat,scale.start=500,shape.start=100))

try(abc(abc2,Nabc=3,probes=probes.good,scale=scale.dat,epsilon=1,
        proposal=function(theta,...)stop("urp!")))
try(abc(pomp(ou2,dprior=function(params,log,...)stop("oof!")),
        Nabc=3,probes=probes.good,scale=scale.dat,epsilon=1,
        proposal=mvn.diag.rw(c(alpha.1=0.01,alpha.2=0.01))))
neval <- 0
try(abc(abc2,Nabc=3,probes=probes.good,scale=scale.dat,epsilon=1,
        proposal=function(theta,...) {
            neval <<- neval+1
            if (neval>1) stop("bif!") else theta
        }))
neval <- 0
try(abc(
    pomp(ou2,dprior=function(params,log,...) {
        neval <<- neval+1
        if (neval>10) stop("eep!")
        else if (log) 0
        else 1
    }),Nabc=20,probes=probes.good,scale=scale.dat,epsilon=1,
    proposal=mvn.diag.rw(c(alpha.1=0.01,alpha.2=0.01))))
neval <- 0
try(abc(ou2,Nabc=10,scale=scale.dat,epsilon=1,
        proposal=mvn.diag.rw(c(alpha.1=0.01,alpha.2=0.01)),
        probes=function(y) {
            neval <<- neval+1
            if (neval>0) stop("pow!")
            else 0
        }))
neval <- 0
try(abc(ou2,Nabc=10,scale=scale.dat,epsilon=1,
        proposal=mvn.diag.rw(c(alpha.1=0.01,alpha.2=0.01)),
        probes=function(y) {
            neval <<- neval+1
            if (neval>1) stop("paf!")
            else 0
        }))

dev.off()

