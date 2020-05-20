\donttest{
sir(
  times=seq(from=0,to=1,by=1/26),
  t0=-1/26
) %>%
  pmcmc(
    Np = 1000,
    Nmcmc = 50,
    dprior = function (beta1, beta2, beta3, ..., log = FALSE) {
      lp <- sum(dnorm(c(beta1,beta2,beta3),mean=400,sd=100,log=TRUE))
      if (log) lp else exp(lp)
    },
    proposal = mvn.diag.rw(c(beta1=1,beta2=1,beta3=1))
  ) -> sirpm

sirpm %>%
  traces(c("beta1","beta2","beta3")) %>%
  plot()
}
