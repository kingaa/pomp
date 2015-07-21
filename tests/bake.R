### test of reproducibility utilities

library(pomp) 

set.seed(5499)
runif(2)
freeze({runif(10)},seed=499586)
runif(2)
freeze(runif(10),seed=499586)
set.seed(5499)
runif(4)

set.seed(5499)
x1 <- bake({runif(4)},file=file.path(tempdir(),"bake1.rds"))
x2 <- bake({runif(4)},file=file.path(tempdir(),"bake2.rds"),seed=5499)
x3 <- bake({runif(4)},file=file.path(tempdir(),"bake1.rds"))
stopifnot(all.equal(x1,x2))
stopifnot(all.equal(x1,x3))

set.seed(5499)
stew({x1 <- runif(4)},file=file.path(tempdir(),"stew1.rds"))
stew({x2 <- runif(4)},file=file.path(tempdir(),"stew2.rds"),seed=5499)
stew({x3 <- runif(4)},file=file.path(tempdir(),"stew1.rds"))
stopifnot(all.equal(x1,x2))
stopifnot(all.equal(x1,x3))
