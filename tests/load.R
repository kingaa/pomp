library(pomp)

pompExample(gompertz)

pomp(gompertz,rmeasure=Csnippet("
  Y = rlnorm(log(X),tau);"),
  statenames="X",paramnames="tau",
  cdir=getwd(),cfile="sf"
) -> po

file.remove(paste0("sf",.Platform$dynlib.ext))

capture.output(x <- simulate(po,verbose=TRUE)) -> out
stopifnot(sum(grepl("loading",out))==2)

pompExample(euler.sir)

solibs(euler.sir) <- NULL
solibs(euler.sir) <- euler.sir@solibs[[1]]
stopifnot(length(euler.sir@solibs)==2,
  length(unique(sapply(euler.sir@solibs,getElement,"src")))==1)
