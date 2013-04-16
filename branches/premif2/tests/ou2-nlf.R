library(pomp)

set.seed(594861940L)

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
          trace=1,
          verbose=TRUE,
          eval.only=TRUE,
          seed=384886L,
          lql.frac = 0.025
          )
se1 <- rep(NA,length(estnames))
print(m1)

m2 <- nlf(
          object=ou2,
          start=theta.truth,
          est=estnames,
          lags=c(4,6),
          nconverge=100,
          nasymp=2000,
          method="Nelder-Mead",
          maxit=500,
          trace=1,
          verbose=TRUE,
          seed=384886L,
          lql.frac = 0.025
          )

se2 <- m2$se

m3 <- nlf(
          object=ou2,
          start=theta.guess,
          lags=c(4,6),
          nconverge=100,
          nasymp=2000,
          method="Nelder-Mead",
          maxit=500,
          trace=1,
          verbose=TRUE,
          eval.only=TRUE,
          seed=384886L,
          lql.frac = 0.025
          )
se3 <- rep(NA,length(estnames))

m4 <- nlf(
          object=ou2,
          start=theta.guess,
          est=estnames,
          lags=c(4,6),
          nconverge=100,
          nasymp=2000,
          method="Nelder-Mead",
          maxit=500,
          trace=1,
          verbose=TRUE,
          seed=384886L,
          lql.frac = 0.025
          )

se4 <- m4$se

print(
      cbind(
            guess=c(theta.guess[estnames],se=se3,value=m3),
            truth=c(theta.truth[estnames],se=se1,value=m1),
            fit.from.guess=c(m4$params[estnames],se=se4,value=m4$value),
            fit.from.truth=c(m2$params[estnames],se=se2,value=m2$value)
            )
      )
