library(pomp.devel)

data(ou2)

set.seed(1066)

po <- ou2
coef(po,c("x1.0","x2.0","alpha.1","alpha.4")) <- c(0,0,0.1,0.2)
po <- simulate(po,times=(1:10000))[[1]]
p.truth <- coef(po)
guess <- p.truth
## guess[c('x1.0','x2.0','alpha.1','alpha.4')] <- c(45,-60,0.8,0.9)

m1 <- nlf(
          object=po,
          start=guess,
          est=c("alpha.1","alpha.4"),
          lags=c(1,2),
          nconverge=100,
          nasymp=2500,
          method="Nelder-Mead",
          maxit=500,
          trace=1,
          verbose=TRUE,
          lql.frac = 0.025
          )
