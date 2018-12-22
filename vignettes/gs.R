library(foreach)
library(doMPI)
library(doRNG)
cl <- startMPIcluster()
registerDoMPI(cl)
registerDoRNG(348885445L)


source("getting_started.R",echo=TRUE)


closeCluster(cl)
try(detach("package:doMPI",unload=TRUE),silent=TRUE)
if (exists("mpi.exit")) mpi.exit()
try(detach("package:Rmpi",unload=TRUE),silent=TRUE)

