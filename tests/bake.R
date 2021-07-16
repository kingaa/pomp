### test of reproducibility utilities

library(pomp)
options(digits=2)

set.seed(5499)
w1 <- runif(2)
w4 <- freeze({runif(5)},seed=c(499586,588686,39995866))
w2 <- runif(2)
w5 <- freeze(runif(5),seed=499586)
set.seed(5499)
w3 <- runif(4)
stopifnot(
  all.equal(c(w1,w2),w3),
  all(w4==w5)
)

set.seed(32765883)
x1 <- bake({runif(4)},file=file.path(tempdir(),"bake1.rds"))
x2 <- bake({runif(4)},file=file.path(tempdir(),"bake2.rds"),seed=32765883)
x3 <- bake({runif(4)},file=file.path(tempdir(),"bake1.rds"))
x3a <- bake({  runif(4)},file=file.path(tempdir(),"bake1.rds"))
rm(.Random.seed)
x4 <- bake({runif(4)},file=file.path(tempdir(),"bake3.rds"),seed=59566)
x5 <- bake({runif(5)},file=file.path(tempdir(),"bake1.rds"))
x6 <- bake({runif(5)},file=file.path(tempdir(),"bake1.rds"))
x7 <- bake({x1+runif(1)},file=file.path(tempdir(),"bake4.rds"))
x8 <- bake({x1+runif(1)},
  file=file.path(tempdir(),"bake4.rds"),
  dependson=x1)
x9 <- bake({x1+runif(1)},
  file=file.path(tempdir(),"bake4.rds"),
  dependson=x1)
x10 <- bake({x1+runif(1)},
  file=file.path(tempdir(),"bake4.rds"),
  dependson=c(x1,x6))
c <- function(x) x+5
x11 <- bake({x1+runif(1)},
  file=file.path(tempdir(),"bake4.rds"),
  dependson=c(x1,c))
x12 <- bake({x1+runif(1)},
  file=file.path(tempdir(),"bake4.rds"),
  dependson=list(c,x1,x13))
stopifnot(
  all.equal(as.numeric(x1),as.numeric(x2)),
  identical(x1,x3),
  identical(x3,x3a),
  identical(x5,x6),
  !identical(x3,x5),
  !identical(x7,x8),
  identical(x8,x9),
  !identical(x9,x10),
  identical(x11,x12)
)

set.seed(113848)
stew({y1 <- runif(4)},file=file.path(tempdir(),"stew1.rda"))
stew({y2 <- runif(4)},file=file.path(tempdir(),"stew2.rda"),seed=113848)
y3 <- y1
print(stew({y1 <- runif(4)},file=file.path(tempdir(),"stew1.rda")))
rm(.Random.seed)
stew({y4 <- runif(4)},file=file.path(tempdir(),"stew3.rda"),seed=59566)
y5 <- y1+y2 
stew({y6 <- y1+y2},dependson=c(y1,y2),file=file.path(tempdir(),"stew3.rda"))
stopifnot(
  identical(y1,y2),
  identical(y1,y3),
  identical(y5,y6)
)
y1 <- 0
stew({y6 <- y1+y2},dependson=c(y1,y2),file=file.path(tempdir(),"stew3.rda"))
stopifnot(identical(y2,y6))

window(sir2(),end=0.5) -> po
simulate(po,seed=1347484107L) -> x
freeze(simulate(po),seed=1347484107L) -> y
attr(y,"seed") <- NULL
attr(y,"kind") <- NULL
attr(y,"normal.kind") <- NULL
stopifnot(identical(x,y))

stopifnot(
  is.null(freeze({rnorm(5); NULL},seed=3494995)),
  is.character(bake({rnorm(5); NULL},seed=3494995,
    file=file.path(tempdir(),"bake4.rds")))
)

rm(.Random.seed)
invisible(bake({runif(4)},file=file.path(tempdir(),"b99.rds"),seed=32765883))
rm(.Random.seed)
invisible(stew({runif(4)},file=file.path(tempdir(),"s99.rda"),seed=32765883))
rm(.Random.seed)
invisible(freeze({runif(4)},seed=32765883))
invisible(freeze({runif(4)}))
