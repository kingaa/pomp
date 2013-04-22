library(pompExamples)

all <- c("food","para1","para2","tri")

sapply(all,function(n)eval(bquote(budmoth.sim(.(n))))) -> bm

names(bm)
x <- lapply(bm,as,"data.frame")

print(lapply(x,tail))

y <- simulate(budmoth.sim(food),seed=3434996L,as.data.frame=TRUE)
tail(y)

z <- trajectory(budmoth.sim(tri),as.data.frame=TRUE)
tail(z)

pf <- pfilter(budmoth.sim(food),seed=34348885L,Np=1000)
logLik(pf)
