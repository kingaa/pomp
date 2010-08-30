library(pomp)

set.seed(594861940L)

data(ou2)

p.truth <- coef(ou2)

m1 <- nlf(
          object=ou2,
          start=p.truth,
          lags=c(4,6),
          nconverge=100,
          nasymp=2000,
          trace=1,
          verbose=TRUE,
          eval.only=TRUE,
          seed=384886L,
          lql.frac = 0.025
          )
print(m1)

m2 <- nlf(
          object=ou2,
          start=p.truth,
          est=c("alpha.2","alpha.3","tau"),
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

se <- m2$se
print(cbind(truth=p.truth[names(se)],fit=m2$params[names(se)],se=se),digits=3)
