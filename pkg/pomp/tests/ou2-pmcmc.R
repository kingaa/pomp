if (Sys.getenv("POMP_FULL_TESTS")=="yes") {

  library(pomp)

  pompExample(ou2)

  dprior.ou2 <- function (params, log, ...) {
    f <- sum(dunif(params,min=coef(ou2)-1,max=coef(ou2)+1,log=TRUE))
    if (log) f else exp(f)
  }

  pdf(file="ou2-pmcmc.pdf")

  f1 <- pmcmc(
              pomp(ou2,dprior=dprior.ou2),
              Nmcmc=20,
              rw.sd=c(alpha.2=0.001,alpha.3=0.001),
              Np=100,
              max.fail=100, 
              verbose=FALSE
              )
  f1 <- continue(f1,Nmcmc=200,max.fail=100)
  plot(f1)

  ff <- pfilter(f1)
  f2 <- pmcmc(
              ff,
              Nmcmc=20,
              rw.sd=c(alpha.2=0.01,alpha.3=0.01),
              max.fail=100, 
              verbose=FALSE
              )

  f3 <- pmcmc(
              ff,
              Nmcmc=20,
              transform=TRUE,
              rw.sd=c(alpha.2=0.01,alpha.3=0.01),
              max.fail=100, 
              verbose=FALSE
              )
  f4 <- pmcmc(f3)
  f4 <- continue(f4,Nmcmc=100)

  if (FALSE) {
    f2 <- pmcmc(
                f1,Nmcmc=1000,Np=500,max.fail=100,
                verbose=FALSE
                )
    plot(f2)
    runs <- rle(conv.rec(f2)[,'loglik'])$lengths
    plot(runs)
    acf(conv.rec(f2)[,c("alpha.2","alpha.3")])
  }

  dprior.ou2 <- function (params, log, ...) {
    f <- sum(dnorm(params,mean=coef(ou2),sd=1,log=TRUE))
    if (log) f else exp(f)
  }

  f5 <- pmcmc(
              ou2,
              start=coef(ou2),
              Nmcmc=20,
              rw.sd=c(alpha.2=0.001,alpha.3=0.001),
              Np=100,
              max.fail=100, 
              verbose=FALSE
              )
  f5 <- continue(f5,Nmcmc=200,max.fail=100)
  plot(f5)

  dev.off()

}
