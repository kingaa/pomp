library(pomp)

set.seed(54588699L)

pompExample(ricker)
coef(ricker,"sigma") <- 0
tm <- sort(runif(n=20,max=3))
x <- trajectory(ricker,times=tm)["N",,]
y <- simulate(ricker,times=tm,states=TRUE)["N",,]
stopifnot(identical(x,y))

pompExample(euler.sir)
tm <- sort(runif(n=100,max=1))
x <- trajectory(euler.sir,times=tm)["I",,]
y <- simulate(euler.sir,times=tm,states=TRUE)["I",,]
table(cut(x-y,breaks=c(-Inf,seq(-0.2,0.2,by=0.01),Inf),ordered=T))
