library(pomp)

pompExample("ou2")

dprior.ou2 <- function (params, hyperparams, ..., log) {
  f <- sum(dunif(params,min=hyperparams$min,max=hyperparams$max,log=TRUE))
  if (log) f else exp(f)
}

pdf(file="ou2-pmcmc.pdf")

f1 <- pmcmc(
            ou2,
            start=coef(ou2),
            Nmcmc=20,
            dprior=dprior.ou2,
            hyperparams=list(min=coef(ou2)-1,max=coef(ou2)+1),
            rw.sd=c(alpha.2=0.01,alpha.3=0.01),
            Np=100,
            max.fail=100, 
            verbose=FALSE
            )
f1 <- continue(f1,Nmcmc=20,max.fail=100)
plot(f1)

ff <- pfilter(f1)
f2 <- pmcmc(
            ff,
            Nmcmc=20,
            dprior=dprior.ou2,
            hyperparams=list(min=coef(ou2)-1,max=coef(ou2)+1),
            rw.sd=c(alpha.2=0.01,alpha.3=0.01),
            max.fail=100, 
            verbose=FALSE
            )


f3 <- pmcmc(
            ff,
            Nmcmc=20,
            dprior=dprior.ou2,
            transform=TRUE,
            hyperparams=list(min=coef(ou2)-1,max=coef(ou2)+1),
            rw.sd=c(alpha.2=0.01,alpha.3=0.01),
            max.fail=100, 
            verbose=FALSE
            )


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

dev.off()
