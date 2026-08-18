library(pomp)

simulate(
  t0=0,
  times=1:10,
  statenames="x",
  obsnames="y",
  params=c(x_0=0),
  rprocess=onestep(function(...)c(x=1)),
  rmeasure=Csnippet("")
) -> s1
x <- states(s1)
y <- obs(s1)
stopifnot(
  all(x==1),
  all(is.na(y)),
  dim(x)==c(1,10),
  dim(y)==c(1,10)
)

simulate(
  t0=0,
  times=1:10,
  statenames="x",
  params=c(x_0=0),
  rprocess=onestep(function(...)c(x=1))
) -> s2
x <- states(s2)
y <- obs(s2)
stopifnot(
  all(x==1),
  dim(x)==c(1,10),
  dim(y)==c(0,10)
)
