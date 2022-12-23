\donttest{

  if (require(subplex)) {

    ricker() |>
      nlf_objfun(est=c("r","sigma","N_0"),lags=c(4,6),
        partrans=parameter_trans(log=c("r","sigma","N_0")),
        paramnames=c("r","sigma","N_0"),
        ti=100,tf=2000,seed=426094906L) -> m1

    subplex(par=log(c(20,0.5,5)),fn=m1,control=list(reltol=1e-4)) -> out

    m1(out$par)
    coef(m1)
    plot(simulate(m1))

  }
}
