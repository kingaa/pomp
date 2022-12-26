library(pomp)

ou2() -> po

try(
  dprocess(po,x=states(po)[,c(1:9,15)],times=time(po)[c(1:9,15)])
)
