if (require(tidyr) && require(ggplot2)) {

  bsflu |>
    gather(variable,value,-date,-day) |>
    ggplot(aes(x=date,y=value,color=variable))+
    geom_line()+
    labs(y="number of boys",title="boarding school flu outbreak")+
    theme_bw()

}
