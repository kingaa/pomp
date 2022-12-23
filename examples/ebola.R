if (require(ggplot2) && require(tidyr)) {

  data(ebolaWA2014)

  ebolaWA2014 |>
    gather(variable,count,cases,deaths) |>
    ggplot(aes(x=date,y=count,group=country,color=country))+
    geom_line()+
    facet_grid(variable~.,scales="free_y")+
    theme_bw()+
    theme(axis.text=element_text(angle=-90))

  ebolaWA2014 |>
    gather(variable,count,cases,deaths) |>
    ggplot(aes(x=date,y=count,group=variable,color=variable))+
    geom_line()+
    facet_grid(country~.,scales="free_y")+
    theme_bw()+
    theme(axis.text=element_text(angle=-90))

  plot(ebolaModel(country="SLE"))
  plot(ebolaModel(country="LBR"))
  plot(ebolaModel(country="GIN"))
}
