library(pomp)

set.seed(1420306530L)

data(dacca)

x <- as.data.frame(dacca)
print(names(x))
print(dim(x))

x <- simulate(dacca,nsim=3,as.data.frame=TRUE)
print(names(x))
print(dim(x))
