\donttest{

  ricker() %>%
    traj_objfun(
      est=c("r","sigma","N_0"),
      partrans=parameter_trans(log=c("r","sigma","N_0")),
      paramnames=c("r","sigma","N_0"),
      ) -> f

  f(log(c(20,0.3,10)))

  if (require(subplex)) {
    subplex(fn=f,par=log(c(20,0.3,10)),control=list(reltol=1e-5)) -> out
  } else {
    optim(fn=f,par=log(c(20,0.3,10)),control=list(reltol=1e-5)) -> out
  }

  f(out$par)

  if (require(ggplot2)) {

    f %>%
      trajectory(format="data.frame") %>%
      ggplot(aes(x=time,y=N))+geom_line()+theme_bw()

  }

}
