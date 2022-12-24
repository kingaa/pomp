options(digits=3,dplyr.summarise.inform=FALSE)
png(filename="pfilter2-%02d.png",res=100)

library(pomp)
library(dplyr)
library(tidyr)
library(ggplot2)

set.seed(9994847L)

ou2(tau=5) %>%
  window(end=10) %>%
  pfilter(Np=5000,save.states="weighted") -> pf

pf %>% saved_states(format="d") %>% names()
pf %>% saved_states(format="d") %>% dim()

pf %>%
  saved_states(format="d") %>%
  pivot_wider(names_from="variable") %>%
  group_by(time) %>%
  summarize(
    p=c(0.05,0.5,0.95),
    x1=wquant(x1,weights=exp(.log.weight),probs=p),
    x2=wquant(x2,weights=exp(.log.weight),probs=p)
  ) %>%
  ungroup() %>%
  pivot_longer(c(x1,x2)) %>%
  pivot_wider(names_from=p) %>%
  ggplot(aes(x=time,y=`0.5`,ymin=`0.05`,ymax=`0.95`,group=name))+
  geom_ribbon(alpha=0.5)+
  geom_line(color='red')+
  facet_grid(name~.)+
  labs(y="")
warnings()

dev.off()
