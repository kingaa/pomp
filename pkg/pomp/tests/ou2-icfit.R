library(pomp)

set.seed(921625222L)

pompExample("ou2")

ics <- c("x1.0","x2.0")    # names of the initial condition parameters

p <- coef(ou2)
p[ics] <- c(4,-4)

fit <- mif(
           ou2,
           start=p,
           Nmif=2,
           ic.lag=1000,
           var.factor=4,
           ivps=ics,
           rw.sd=c(
             x1.0=1,x2.0=1
             ),
           Np=1000,
           cooling.factor=1,
           max.fail=10
           )

fit <- mif(
           ou2,
           start=p,
           Nmif=2,
           ic.lag=10,
           var.factor=4,
           ivps=ics,
           rw.sd=c(
             x1.0=1,x2.0=1
             ),
           Np=1000,
           cooling.factor=1,
           max.fail=10
           )

fit <- mif(
           window(ou2,end=10),
           start=p,
           Nmif=500,
           ic.lag=10,
           var.factor=4,
           ivps=ics,
           rw.sd=c(
             x1.0=1,x2.0=1
             ),
           Np=1000,
           cooling.factor=1,
           max.fail=10
           )

print(logLik(pfilter(ou2,Np=2000)),digits=4)
print(logLik(pfilter(ou2,params=p,Np=2000)),digits=4)
print(logLik(pfilter(ou2,params=coef(fit),Np=2000)),digits=4)
print(coef(fit,ics))
print(coef(ou2,ics))
print(p-coef(ou2))
print(coef(fit)-p)
print(coef(fit)-coef(ou2))
