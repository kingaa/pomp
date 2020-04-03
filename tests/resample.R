library(pomp)
set.seed(807969746L)

w <- runif(100)
k <- systematic_resample(w,100)
try(k <-systematic_resample(-w))
k <- systematic_resample(w,10)
k <- systematic_resample(w,110)
