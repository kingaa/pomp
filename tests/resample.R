library(pomp)
library(magrittr)
set.seed(807969746L)

w <- runif(100)
k <- .Call(pomp:::P_systematic_resampling,w)
try(k <- .Call(pomp:::P_systematic_resampling,-w))
