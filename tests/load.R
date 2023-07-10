library(pomp)

gompertz() |>
  pomp(
    rmeasure=Csnippet(r"{
  Y = rlnorm(log(X),Tau);
  }"),
  statenames="X",paramnames="tau",
  globals=Csnippet("static double Tau = 0;"),
  on_load=Csnippet("Tau = 4.3;"),
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
