library(pomp)

ou2() -> po

stopifnot(
  round(sum(dprocess(po)),3)==1.400
)

try(
  dprocess(po,x=states(po)[,c(1:9,15)],times=time(po)[c(1:9,15)])
)
