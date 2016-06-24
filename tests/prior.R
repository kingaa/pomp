library(pomp)

pompExample(ou2)
coef(ou2,"alpha.sd") <- 5

set.seed(1835425749L)

po <- pomp(ou2,
           dprior=function(params,log,...) {
             ll <- sum(
                       dnorm(
                             x=params[c("alpha.1","alpha.2","alpha.3","alpha.4")],
                             mean=c(0.8,-0.5,0.3,0.9),
                             sd=params["alpha.sd"],
                             log=TRUE
                             )
                       )
             if (log) ll else exp(ll)
           },
           rprior=function(params,...) {
             params[c("alpha.1","alpha.2","alpha.3","alpha.4")] <- rnorm(
                                                                         n=4,
                                                                         mean=c(0.8,-0.5,0.3,0.9),
                                                                         sd=params["alpha.sd"]
                                                                         )
             params
           }
           )


stopifnot(all.equal(mean(dprior(po,params=parmat(coef(po),3))),dnorm(x=0,mean=0,sd=5)^4))
rprior(po,params=coef(po))

coef(po,"alpha.sd") <- 1
mean(dprior(po,params=rprior(po,params=parmat(coef(po),10000)),log=TRUE))+0.5*(1+log(2*pi))*4

po <- pomp(po,
           rprior=Csnippet("
              alpha_1 = rnorm(0.8,alpha_sd);
              alpha_2 = rnorm(-0.5,alpha_sd);
              alpha_3 = rnorm(0.3,alpha_sd);
              alpha_4 = rnorm(0.9,alpha_sd);"),
           dprior=Csnippet("
              lik = 0;
              lik += dnorm(alpha_1,0.8,alpha_sd,1);
              lik += dnorm(alpha_2,-0.5,alpha_sd,1);
              lik += dnorm(alpha_3,0.3,alpha_sd,1);
              lik += dnorm(alpha_4,0.9,alpha_sd,1);
              lik = (give_log) ? lik : exp(lik);"),
           paramnames=c("alpha.1","alpha.2","alpha.3","alpha.4","alpha.sd"))
rprior(po,params=coef(po))
mean(dprior(po,params=rprior(po,params=parmat(coef(po),10000))))

po <- pomp(po,rprior=function(params,...)runif(length(params)))
rprior(po,params=coef(po))

po <- pomp(po,rprior=function(params,...)runif(2))
try(rprior(po,params=coef(po)))
