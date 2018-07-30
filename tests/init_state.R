library(pomp)
library(magrittr)
set.seed(807969746L)

pompExample(gompertz)
init.state(gompertz)
init.state(gompertz,params=coef(gompertz))

p <- coef(gompertz)[-5]
try(init.state(gompertz,params=p))

gompertz %>% simulate(initializer=NULL)

gompertz %>%
  pomp(initializer=function (params, t0, ...) 5) -> po
try(init.state(po))

pompExample(sir)
try(sir %>% simulate(initializer=NULL))
