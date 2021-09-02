\dontrun{
  parus %>%
    pfilter(Np=1000,times="year",t0=1960,
      params=c(K=190,r=2.7,sigma=0.2,theta=0.05,N.0=148),
      rprocess=discrete_time(
        function (r, K, sigma, N, ...) {
          e <- rnorm(n=1,mean=0,sd=sigma)
          c(N = exp(log(N)+r*(1-N/K)+e))
        },
        delta.t=1
      ),
      rmeasure=function (N, theta, ...) {
        c(pop=rnbinom(n=1,size=1/theta,mu=N+1e-10))
      },
      dmeasure=function (pop, N, theta, ..., log) {
        dnbinom(x=pop,mu=N+1e-10,size=1/theta,log=log)
      },
      partrans=parameter_trans(log=c("sigma","theta","N_0","r","K")),
      paramnames=c("sigma","theta","N_0","r","K")
    ) -> pf

  pf %>% logLik()

  pf %>% simulate() %>% plot()
}
