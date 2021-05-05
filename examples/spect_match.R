\donttest{

  ricker() %>%
    spect_objfun(
      est=c("r","sigma","N_0"),
      partrans=parameter_trans(log=c("r","sigma","N_0")),
      paramnames=c("r","sigma","N_0"),
      kernel.width=3,
      nsim=100,
      seed=5069977
    ) -> f

  f(log(c(20,0.3,10)))
  f %>% spect() %>% plot()

  library(subplex)
  subplex(fn=f,par=log(c(20,0.3,10)),control=list(reltol=1e-5)) -> out
  f(out$par)

  f %>% summary()

  f %>% spect() %>% plot()

}
