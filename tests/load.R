library(pomp)

gompertz() -> gompertz

pomp(gompertz,rmeasure=Csnippet("
  Y = rlnorm(log(X),tau);"),
  statenames="X",paramnames="tau",
  cdir=getwd(),cfile="sf"
) -> po

file.remove(paste0("sf",.Platform$dynlib.ext))

capture.output(x <- simulate(po,verbose=TRUE)) -> out
stopifnot(sum(grepl("loading",out))==2)

sir() -> sir

solibs(sir) <- NULL
solibs(sir) <- sir@solibs[[1]]
stopifnot(length(sir@solibs)==2,
  length(unique(sapply(sir@solibs,getElement,"src")))==1)
