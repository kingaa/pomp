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
            Nabc=10000,
            start=coef(ou2),
            pars=c("alpha.1","alpha.2"),
            probes=probes.good,
            scale=scale.dat,
            epsilon=1.7,
            rw.sd= c(alpha.1=0.01,alpha.2=0.01)
            )

plot(abc1,scatter=TRUE)
plot(abc1)

## check how sticky the chain is:
runs <- rle(as.vector(conv.rec(abc1)[, "alpha.1"]))
hist(runs$lengths)
mean(runs$length)

abc2 <- abc(po,
            Nabc=2000,
            pars=c("alpha.1","alpha.2"),
            probes=probes.good,
            scale=scale.dat,
            epsilon=1,
            rw.sd= c(alpha.1=0.01,alpha.2=0.01)
            )
plot(abc2)

abc3 <- abc(po,
            Nabc=2000,
            probes=probes.good,
            scale=scale.dat,
            epsilon=2,
            rw.sd= c(alpha.1=0.01,alpha.2=0.01)
            )
abc3 <- continue(abc3,Nabc=3000)
plot(abc3)

abc4 <- abc(probe(po,probes=probes.good,nsim=200),
            Nabc=2000,
            scale=scale.dat,
            epsilon=2,
            rw.sd= c(alpha.1=0.01,alpha.2=0.01)
            )
plot(abc4)

abc5 <- abc(abc4,Nabc=1000)
plot(abc5)

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
            Nabc=2000,
            pars=c("alpha.1","alpha.2"),
            probes=probes.good,
            scale=scale.dat,
            epsilon=1,
            rw.sd= c(alpha.1=0.01,alpha.2=0.01)
            )
plot(abc6)

try(abc7 <- c(abc2,abc3))
plot(abc7 <- c(abc2,abc4))
plot(abc7,scatter=TRUE)
plot(conv.rec(c(abc2,abc4)))
plot(conv.rec(c(abc7,abc6)))

dev.off()

