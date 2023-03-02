\dontrun{
if (require(ggplot2) && require(tidyr)) {

  ebolaWA2014 |>
    pivot_longer(c(cases,deaths)) |>
    ggplot(aes(x=date,y=value,group=name,color=name))+
    geom_line()+
    facet_grid(country~.,scales="free_y")+
    theme_bw()+
    theme(axis.text=element_text(angle=-90))

}

plot(ebolaModel(country="SLE"))
plot(ebolaModel(country="GIN"))
plot(ebolaModel(country="LBR"))
}
