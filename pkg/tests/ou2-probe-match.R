library(pomp)
set.seed(1066L)

data(ou2)
po <- ou2
coef(po,c("x1.0","x2.0","alpha.1","alpha.4")) <- c(0,0,0.1,0.2)

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
               po,
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
                  y1.lag3=probe.acf(var="y1",lag=3,type="corr"),
                  y2.cov12=probe.cov(vars=c("y1"),lag=12,method="spearman"),
                  y12.cov8=probe.cov(vars=c("y2","y1"),lag=8,method="pearson")
                  ),
                nsim=500
                )
                
summary(pm.ou2)
