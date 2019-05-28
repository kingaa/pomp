library(pomp)
set.seed(807969746L)

gompertz() -> gompertz
rinit(gompertz)
rinit(gompertz,params=coef(gompertz))

p <- coef(gompertz)[-5]
try(rinit(gompertz,params=p))

gompertz %>% simulate(rinit=NULL)

gompertz %>%
  pomp(rinit=function (...) 5) -> po
try(rinit(po))
try(rinit("gompertz"))
try(rinit())

pp <- parmat(coef(gompertz),10)
stopifnot(gompertz %>% rinit(params=pp) %>% as.numeric()==1)
rinit(gompertz,params=pp[,1:3],nsim=2) -> x0
stopifnot(dim(x0)==c(1,6))
dimnames(x0)
colnames(pp) <- head(LETTERS,10)
rinit(gompertz,params=pp[,1:5],nsim=2) -> x0
stopifnot(dim(x0)==c(1,10))
dimnames(x0)
rinit(gompertz,params=pp[,1:5],nsim=1) -> x0
stopifnot(dim(x0)==c(1,5),colnames(x0)==head(LETTERS,5))

try(gompertz %>%
  pomp(rinit=function(...)
    c(r=32)) %>%
  rinit())
try({
  pp <- matrix(c(1:5),1,5)
  rownames(pp) <- "a"
  gompertz %>%
    pomp(rinit=function(a,...)
      c(X=rep(1,a))) %>%
    rinit(params=pp)
})

sir() -> sir
try(sir %>% simulate(rinit=NULL))
sir %>%
  pomp(rinit=function(seas_1,...)
    c(S=seas_1)) %>%
  rinit()

gompertz %>% rinit(nsim=3) -> x
gompertz %>% pomp(rinit=function(K,...)c(X=K)) %>% rinit(nsim=3) -> y
stopifnot(identical(x,y))
