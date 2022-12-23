\donttest{

  if (require(ggplot2)) {

    rw2() |> plot()

    rw2(s1=1,s2=1,tau=0.1) |>
      simulate(nsim=10,format="d") |>
      ggplot(aes(x=y1,y=y2,group=.id,color=.id))+
      geom_path()+
      guides(color="none")+
      theme_bw()

  }
}
