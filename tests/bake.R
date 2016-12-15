### test of reproducibility utilities

library(pomp)

set.seed(5499)
w1 <- runif(2)
freeze({runif(5)},seed=499586)
w2 <- runif(2)
freeze(runif(5),seed=499586)
set.seed(5499)
w3 <- runif(4)
stopifnot(all.equal(c(w1,w2),w3))

set.seed(32765883)
x1 <- bake({runif(4)},file=file.path(tempdir(),"bake1.rds"))
x2 <- bake({runif(4)},file=file.path(tempdir(),"bake2.rds"),seed=32765883)
x3 <- bake({runif(4)},file=file.path(tempdir(),"bake1.rds"))
stopifnot(all.equal(as.numeric(x1),as.numeric(x2)))
stopifnot(all.equal(as.numeric(x1),as.numeric(x3)))

set.seed(113848)
stew({y1 <- runif(4)},file=file.path(tempdir(),"stew1.rds"))
stew({y2 <- runif(4)},file=file.path(tempdir(),"stew2.rds"),seed=113848)
print(stew({y3 <- runif(4)},file=file.path(tempdir(),"stew1.rds")))
stopifnot(all.equal(y1,y2))
try(stopifnot(all.equal(y1,y3)))

pompExample(gillespie.sir)
simulate(gillespie.sir,seed=1347484107L) -> x
freeze(simulate(gillespie.sir),seed=1347484107L) -> y
stopifnot(identical(x,y))

detach("package:pomp", unload=TRUE)
