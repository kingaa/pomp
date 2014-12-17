library(pompExamples)

all <- c("food","para1","para2","tri")

bm <- pompExample(budmoth,envir=NULL)

names(bm)
x <- lapply(bm,as,"data.frame")

print(lapply(x,tail))

y <- simulate(bm$food,seed=3434996L,as.data.frame=TRUE)
tail(y)

z <- trajectory(bm$tri,as.data.frame=TRUE)
tail(z)

pf <- pfilter(bm$para1,seed=34348885L,Np=1000)
logLik(pf)
