\donttest{
  library(magrittr)

  gompertz() -> po
  
  ## A list of probes:
  plist <- list(
    mean=probe.mean("Y",trim=0.1,transform=sqrt),
    sd=probe.sd("Y",transform=sqrt),
    probe.marginal("Y",ref=obs(po)),
    probe.acf("Y",lags=c(1,3,5),type="correlation",transform=sqrt),
    probe.quantile("Y",prob=c(0.25,0.75),na.rm=TRUE)
  )

  ## Construct the probe-matching objective function.
  ## Here, we just want to estimate 'K'.
  po %>%
    probe.objfun(probes=plist,nsim=100,seed=5069977,
      est="K") -> f

  ## Any numerical optimizer can be used to minimize 'f'.
  library(subplex)

  subplex(fn=f,par=0.4,control=list(reltol=1e-5)) -> out

  ## Call the objective one last time on the optimal parameters:
  f(out$par)

  ## There are 'plot' and 'summary' methods:
  f %>% as("probed_pomp") %>% plot()
  f %>% summary()

  f %>% probe() %>% plot()

  ## One can modify the objective function with another call
  ## to 'probe.objfun':

  f %>% probe.objfun(est=c("r","K")) -> f1
}
