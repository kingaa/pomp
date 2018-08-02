library(pomp)
library(magrittr)
set.seed(807969746L)

pompExample(gompertz)
rinit(gompertz)
rinit(gompertz,params=coef(gompertz))

p <- coef(gompertz)[-5]
try(rinit(gompertz,params=p))

gompertz %>% simulate(rinit=NULL)

gompertz %>%
  pomp(rinit=function (params, t0, ...) 5) -> po
try(rinit(po))

pp <- parmat(coef(gompertz),10)
stopifnot(gompertz %>% rinit(params=pp) %>% as.numeric()==1)
try(rinit(gompertz,params=pp,nsim=5))
try(gompertz %>%
  pomp(rinit=function(t0,params,...)
    c(r=32)) %>%
  rinit())
try({
  pp <- matrix(c(1:5),1,5)
  rownames(pp) <- "a"
  gompertz %>%
    pomp(rinit=function(t0,params,...)
      c(X=rep(1,params["a"]))) %>%
    rinit(params=pp)
})


pompExample(sir)
try(sir %>% simulate(rinit=NULL))
sir %>%
  pomp(rinit=function(params,t0,covars,...)
    c(S=covars["seas1"])) %>%
  rinit()
