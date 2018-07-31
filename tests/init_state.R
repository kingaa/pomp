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

pp <- parmat(coef(gompertz),10)
stopifnot(gompertz %>% init.state(params=pp) %>% as.numeric()==1)
try(init.state(gompertz,params=pp,nsim=5))
try(gompertz %>%
  pomp(initializer=function(t0,params,...)
    c(r=32)) %>%
  init.state())
try({
  pp <- matrix(c(1:5),1,5)
  rownames(pp) <- "a"
  gompertz %>%
    pomp(initializer=function(t0,params,...)
      c(X=rep(1,params["a"]))) %>%
    init.state(params=pp)
})


pompExample(sir)
try(sir %>% simulate(initializer=NULL))
sir %>%
  pomp(initializer=function(params,t0,covars,...)
    c(S=covars["seas1"])) %>%
  init.state()
