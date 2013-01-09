library(pomp)

set.seed(54588699L)

pompExample("ricker")
coef(ricker,"sigma") <- 0
tm <- sort(runif(n=20,max=3))
x <- trajectory(ricker,times=tm)["N",,]
y <- simulate(ricker,times=tm,states=TRUE)["N",,]
stopifnot(identical(x,y))

pompExample("verhulst")
coef(verhulst,c("n.0","sigma")) <- c(15,0)
tm <- sort(runif(n=100,max=1))
x <- trajectory(verhulst,times=tm)["n",,]
y <- simulate(verhulst,times=tm,states=TRUE)["n",,]
table(cut(x-y,breaks=10))
