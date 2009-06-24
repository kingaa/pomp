library(pomp)

data(ou2)

set.seed(1066L)

po <- ou2
coef(po,c("x1.0","x2.0","alpha.1","alpha.4")) <- c(0,0,0.1,0.2)
po <- simulate(po)
guess <- p.truth <- coef(po)

m1 <- nlf(
          object=po,
          start=guess,
          est=c("alpha.1","alpha.4"),
          lags=c(1,2),
          nconverge=100,
          nasymp=1000,
          method="Nelder-Mead",
          maxit=500,
          trace=1,
          verbose=TRUE,
          lql.frac = 0.025
          )

se <- p.truth
se[] <- NA
se[names(m1$se)] <- m1$se
print(cbind(truth=p.truth,fit=m1$params,se=se),digits=3)

po <- simulate(po,times=(1:1000))
m2 <- nlf(
          object=po,
          start=guess,
          est=c("alpha.1","alpha.4"),
          lags=c(1,2),
          nconverge=100,
          nasymp=10000,
          method="Nelder-Mead",
          maxit=500,
          trace=1,
          verbose=TRUE,
          lql.frac = 0.025
          )

se <- p.truth
se[] <- NA
se[names(m2$se)] <- m2$se
print(cbind(truth=p.truth,fit=m2$params,se=se),digits=3)
