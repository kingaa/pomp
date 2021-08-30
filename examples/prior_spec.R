\dontrun{
  ## Starting with an existing pomp object

  verhulst() %>% window(end=30) -> po
  
  ## we add or change prior distributions using the two
  ## arguments 'rprior' and 'dprior'. Here, we introduce
  ## a Gamma prior on the 'r' parameter.
  ## We construct 'rprior' and 'dprior' using R functions.

  po %>%
    bsmc2(
      rprior=function (n_0, K0, K1, sigma, tau, r0, r1, ...) {
        c(
          n_0 = n_0,
          K = rgamma(n=1,shape=K0,scale=K1),
          r = rgamma(n=1,shape=r0,scale=r1),
          sigma = sigma,
          tau = tau
        )
      },
      dprior=function(K, K0, K1, r, r0, r1, ..., log) {
        p <- dgamma(x=c(K,r),shape=c(K0,r0),scale=c(K1,r1),log=log)
        if (log) sum(p) else prod(p)
      },
      params=c(n_0=10000,K=10000,K0=10,K1=1000,
        r=0.9,r0=0.9,r1=1,sigma=0.5,tau=0.3),
      Np=1000
    ) -> B

  ## We can also pass them as C snippets:

  po %>%
    bsmc2(
      rprior=Csnippet("
         K = rgamma(K0,K1);
         r = rgamma(r0,r1);"
      ),
      dprior=Csnippet("
         double lik1 = dgamma(K,K0,K1,give_log);
         double lik2 = dgamma(r,r0,r1,give_log);
         lik = (give_log) ? lik1+lik2 : lik1*lik2;"
      ),
      paramnames=c("K","K0","K1","r","r0","r1"),
      params=c(n_0=10000,K=10000,K0=10,K1=1000,
        r=0.9,r0=0.9,r1=1,sigma=0.5,tau=0.3),
      Np=10000
    ) -> B

  ## The prior is plotted in grey; the posterior, in blue.
  plot(B)

  B %>%
    pmcmc(Nmcmc=100,Np=1000,proposal=mvn.diag.rw(c(r=0.01,K=10))) -> Bb

  plot(Bb,pars=c("loglik","log.prior","r","K"))

}
