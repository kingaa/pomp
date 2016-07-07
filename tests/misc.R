library(pomp)
library(magrittr)
set.seed(807969746L)

w <- runif(100)
k <- .Call(pomp:::systematic_resampling,w)
try(k <- .Call(pomp:::systematic_resampling,-w))

pompExample(euler.sir)
euler.sir %<>% pomp(initializer=NULL)
try(simulate(euler.sir))

pompExample(gompertz)
coef(gompertz) <- coef(gompertz)[-5]
try(simulate(gompertz))
