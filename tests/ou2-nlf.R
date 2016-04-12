library(pomp)

set.seed(583615606L)

pompExample(ou2)
estnames=c("alpha.2","alpha.3","tau")
theta.truth <- coef(ou2)
theta.guess <- theta.truth
theta.guess[estnames] <- theta.guess[estnames]*1.5

m1 <- nlf(
          object=ou2,
          start=theta.truth,
          lags=c(4,6),
          nconverge=100,
          nasymp=2000,
          eval.only=TRUE,
          seed=426094906L,
          lql.frac = 0.025
          )

m2 <- nlf(
          m1,
          est=estnames,
          maxit=500,
          method="Nelder-Mead"
          )

m3 <- nlf(
          object=ou2,
          start=theta.guess,
          lags=c(4,6),
          nconverge=100,
          nasymp=2000,
          maxit=500,
          method="Nelder-Mead",
          eval.only=TRUE,
          seed=426094906L,
          lql.frac = 0.025
          )

m4 <- nlf(
          m3,
          est=estnames,
          seed=300983678L
          )

options(scipen=3)
print(
      signif(
             rbind(
                   fit.from.guess=c(coef(m4,estnames),se=m4$se,value=logLik(m4)),
                   fit.from.truth=c(coef(m2,estnames),se=m2$se,value=logLik(m4)),
                   guess=c(theta.guess[estnames],se=rep(NA,length(estnames)),value=logLik(m3)),
                   truth=c(theta.truth[estnames],se=rep(NA,length(estnames)),value=logLik(m1))
                   ),
             4
             )
      )
