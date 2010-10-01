library(pomp)

data(ricker)

pdf(file="ricker-probe.pdf")

set.seed(6457673L)
z <- as.numeric(data.array(ricker))

po <- ricker

pb <- probe(
            po,
            probes=probe.median("y"),
            nsim=1000,
            seed=838775L
            )
plot(pb)
summary(pb)

pb <- probe(
            po,
            probes=probe.nlar("y",lags=c(1,2,3),powers=c(1,1,1),transform="sqrt"),
            nsim=1000,
            seed=838775L
            )
plot(pb)
summary(pb)

pb <- probe(
            po,
            probes=probe.nlar("y",lags=c(1,2,3),powers=1,transform="sqrt"),
            nsim=1000,
            seed=838775L
            )
plot(pb)
summary(pb)

pb <- probe(
            po,
            probes=probe.nlar("y",lags=1,powers=c(1,2,3),transform="sqrt"),
            nsim=1000,
            seed=838775L
            )
plot(pb)
summary(pb)

pb <- probe(
            po,
            probes=probe.marginal(
              var="y",
              transform=sqrt,
              ref=z,
              diff=1,
              order=3
              ),
            nsim=1000,
            seed=838775L
            )
pb@datvals
summary(pb)
plot(pb)

pb <- probe(
            po,
            probes=list(
              probe.marginal(
                               var="y",
                               transform=sqrt,
                               ref=z,
                               diff=1,
                               order=3
                               ),
              probe.acf(
                        var="y",
                        lag.max=5
                        ),
              mean=probe.mean(var="y",transform=sqrt)
              ),
            nsim=1000,
            seed=838775L
            )
pb@datvals
summary(pb)
plot(pb)

coef(po) <- c(log.r=log(10),log.sigma=log(0.3),log.phi=log(20),N.0=5,e.0=0)

pb <- probe(
            po,
            probes=probe.marginal(
              var="y",
              transform=sqrt,
              ref=z,
              diff=1,
              order=3
              ),
            nsim=1000,
            seed=838775L
            )
pb@datvals
summary(pb)
plot(pb)

system.time(
            pm <- probe.match(
                              pb,
                              est=c("log.r","log.phi","N.0"),
                              parscale=c(0.1,0.1,0.1),
                              nsim=1000,
                              seed=838775L,
                              method="Nelder-Mead",
                              reltol=1e-7,
                              fail.value=1e9
                              )
            )
plot(pm)

cbind(truth=coef(ricker),est=coef(pm),guess=coef(po))

pb <- probe(
            po,
            probes=probe.nlar(
              var="y",
              transform=sqrt,
              lags=1,
              powers=c(1,2,3)
              ),
            nsim=1000,
            seed=838775L
            )
pb@datvals
summary(pb)
plot(pb)

system.time(
            pm <- probe.match(
                              pb,
                              est=c("log.r","log.phi","N.0"),
                              parscale=c(0.1,0.1,0.1),
                              nsim=1000,
                              seed=838775L,
                              method="Nelder-Mead",
                              reltol=1e-7,
                              fail.value=1e9
                              )
            )
plot(pm)

cbind(truth=coef(ricker),est=coef(pm),guess=coef(po))

pb <- probe(
            po,
            probes=probe.marginal(
              var="y",
              transform=sqrt,
              ref=runif(length(time(ricker))),
              diff=2,
              order=3
              ),
            nsim=1000,
            seed=838775L
            )
pb@datvals
summary(pb)
plot(pb)

pb <- probe(
            ricker,
            probes=probe.acf(
              var="y",
              transform=sqrt,
              lag.max=5,
              type="cov"
              ),
            nsim=1000,
            seed=838775L
            )
pb@datvals
summary(pb)
plot(pb)

pb <- probe(
            ricker,
            probes=probe.acf(
              var="y",
              transform=sqrt,
              lag.max=5,
              type="cor"
              ),
            nsim=1000,
            seed=838775L
            )
pb@datvals
summary(pb)
plot(pb)

pb <- probe(
            ricker,
            probes=list(
              v=probe.var("y",transform=sqrt),
              probe.acf(
                        var="y",
                        transform=sqrt,
                        lag.max=2,
                        type="cov"
                        ),
              probe.nlar(
                         var="y",
                         transform=sqrt,
                         lags=c(1,2),
                         powers=1
                         )
              ),
            nsim=1000,
            seed=838775L
            )
pb@datvals
summary(pb)
plot(pb)

try(
    probe(
          ricker,
          probes=list(
            mn=probe.mean("y",transform=sqrt,trim=0.1),
            md=probe.median("y",na.rm=FALSE),
            wacko=function(y)y[sample.int(n=length(y),size=sample.int(n=3,size=1))]
            ),
          nsim=100,
          seed=838775L
          )
    )

try(
    probe(
          ricker,
          probes=list(
            mn=probe.mean("y",transform=sqrt,trim=0.1),
            md=function(y)median(as.numeric(y)),
            wacko=function(y) if (y[1]==69) 1 else c(1,2)
            ),
          nsim=100,
          seed=838775L
          )
    )


try(
    probe(
          ricker,
          probes=list(
            mn=probe.mean("y",transform=sqrt,trim=0.1),
            md=function(y)median(as.numeric(y)),
            wacko=function(y) if (y[29]==98) c(1,2) else 1
            ),
          nsim=20,
          seed=838775L
          )
    )

dev.off()
