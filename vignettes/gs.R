if (file.exists("CLUSTER")) {
  scan("CLUSTER",what=integer(0)) -> ncpu
    library(foreach)
library(doMPI)
library(doRNG)
  cl <- startMPIcluster(ncpu,verbose=TRUE,logdir="/tmp")
  registerDoMPI(cl)
registerDoRNG(348885445L)
source("getting_started.R",echo=TRUE)
}