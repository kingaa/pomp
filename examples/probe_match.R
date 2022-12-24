\donttest{
  gompertz() -> po
  
  ## A list of probes:
  plist <- list(
    mean=probe_mean("Y",trim=0.1,transform=sqrt),
    sd=probe_sd("Y",transform=sqrt),
    probe_marginal("Y",ref=obs(po)),
    probe_acf("Y",lags=c(1,3,5),type="correlation",transform=sqrt),
    probe_quantile("Y",prob=c(0.25,0.75),na.rm=TRUE)
  )

  ## Construct the probe-matching objective function.
  ## Here, we just want to estimate 'K'.
  po %>%
    probe_objfun(probes=plist,nsim=100,seed=5069977,
      est="K") -> f

  ## Any numerical optimizer can be used to minimize 'f'.
  if (require(subplex)) {

    subplex(fn=f,par=0.4,control=list(reltol=1e-5)) -> out

  } else {

    optim(fn=f,par=0.4,control=list(reltol=1e-5)) -> out

  }

  ## Call the objective one last time on the optimal parameters:
  f(out$par)
  coef(f)

  ## There are 'plot' and 'summary' methods:
  f %>% as("probed_pomp") %>% plot()
  f %>% summary()

  ## One can convert an objective function to a data frame:
  f %>% as("data.frame") %>% head()
  f %>% as("probed_pomp") %>% as("data.frame") %>% head()

  f %>% probe() %>% plot()

  ## One can modify the objective function with another call
  ## to 'probe_objfun':

  f %>% probe_objfun(est=c("r","K")) -> f1
  optim(fn=f1,par=c(0.3,0.3),control=list(reltol=1e-5)) -> out
  f1(out$par)
  coef(f1)
}
