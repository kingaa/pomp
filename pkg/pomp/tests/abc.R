### OU2 test of abc for pomp;  nov 2013

library(pomp) 
pompExample(ou2)
plot(ou2)

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

scale.dat <- apply(psim@simvals,2,sd)

po <- ou2 

theta <- coef(ou2)

abc1 <- abc(po,
            Nabc=10000,
            start=coef(ou2),
            pars=c("alpha.1","alpha.2"),
            probes=probes.good,
            scale=scale.dat,
            epsilon=3,
            rw.sd= c(alpha.1=0.01,alpha.2=0.01),
            hyperparams=list(junk=0),
            verbose=TRUE
            )

plot(abc1,scatter=TRUE)
plot(abc1)

## check how sticky the chain is:
runs <- rle(as.vector(conv.rec(abc1)[, "alpha.1"]))
hist(runs$lengths)
mean(runs$length)

dev.off()
