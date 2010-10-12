library(pomp)
set.seed(1066L)

data(ou2)

pm.ou2 <- probe(
                ou2,
                probes=list(
                  y1.mean=probe.mean(var="y1"),
                  y2.mean=probe.mean(var="y2"),
                  y1.sd=probe.sd(var="y1"),
                  y2.sd=probe.sd(var="y2")
                  ),
                nsim=500
                )
pm.po <- probe(
               ou2,
               params=c(
                 alpha.1=0.1,alpha.2=-0.5,alpha.3=0.3,alpha.4=0.2,
                 sigma.1=3,sigma.2=-0.5,sigma.3=2,
                 tau=1,
                 x1.0=0,x2.0=0
                 ),
               probes=list(
                 y1.mean=probe.mean(var="y1"),
                 y2.mean=probe.mean(var="y2"),
                 y1.sd=probe.sd(var="y1"),
                 y2.sd=probe.sd(var="y2")
                 ),
               nsim=500
               )

summary(pm.ou2)
summary(pm.po)

plot(pm.ou2)
plot(pm.po)

pm.ou2 <- probe(
                ou2,
                probes=list(
                  y1acf=probe.acf(var="y1",lag.max=2,type="corr"),
                  y2acf=probe.acf(var=c("y2"),lag.max=2),
                  y12ccf=probe.ccf(var=c("y2","y1"),lags=c(3,8))
                  ),
                nsim=500
                )
                
summary(pm.ou2)

pb <- probe(
            ou2,
            probes=list(
              y1=probe.quantile(var="y1",prob=seq(0.1,0.9,by=0.1)),
              probe.acf(var=c("y1","y2"),lag.max=4,transform=identity),
              pd=probe.period(var="y1",kernel.width=3)
              ),
            nsim=200
            )
summary(pb)

po <- ou2
coef(po,c("alpha.2","alpha.3")) <- c(0,0)
coef(po,c("sigma.2","sigma.1","sigma.3")) <- c(0,0.0,0.0)
coef(po,c("tau")) <- c(0.0)
po <- simulate(po)
pb <- probe(
            po,
            probes=list(
              probe.acf(var=c("y1","y2"),lag.max=1,type="cor"),
              probe.nlar("y1",lags=1,powers=1),
              probe.nlar("y2",lags=1,powers=1)
              ),
            nsim=1000,
            seed=1066L
            )
x <- as.data.frame(po)
x <- sweep(x,2,mean(x))
y1 <- head(x$y1,-1)
z1 <- tail(x$y1,-1)
y2 <- head(x$y2,-1)
z2 <- tail(x$y2,-1)
max(abs(pb@datvals-c(mean(y1*z1)/mean(x$y1^2),mean(y2*z2)/mean(x$y2^2),mean(y1*z1)/mean(y1*y1),mean(y2*z2)/mean(y2*y2))))

po <- simulate(ou2)
pb <- probe(
            po,
            probes=list(
              probe.acf(var=c("y1"),lag.max=2,type="cov"),
              probe.ccf(vars=c("y1","y1"),lags=c(0,1,2))
              ),
            nsim=1000,
            seed=1066L
            )
plot(pb)
summary(pb)

pb <- probe(
            po,
            probes=probe.ccf(vars=c("y1","y2"),lags=c(-5,-3,1,4,8)),
            nsim=1000,
            seed=1066L
            )
plot(pb)
summary(pb)

head(as(pb,"data.frame"))


