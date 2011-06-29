## .onAttach <- function (...) {
##   version <- library(help=pomp)$info[[1]]
##   version <- strsplit(version[pmatch("Version",version)]," ")[[1]]
##   version <- version[nchar(version)>0][2]
##   packageStartupMessage("This is pomp version ",version,"\n")
## }

setGeneric("print",function(x,...)standardGeneric("print"))
setGeneric("plot",function(x,y,...)standardGeneric("plot"))
setGeneric("summary",function(object,...)standardGeneric("summary"))
setGeneric("simulate",function(object,nsim=1,seed=NULL,...)standardGeneric("simulate"))
setGeneric("time",function(x,...)standardGeneric("time"))
setGeneric("coef",function(object,...)standardGeneric("coef"))
setGeneric("logLik",function(object,...)standardGeneric("logLik"))
setGeneric("window",function(x,...)standardGeneric("window"))
setGeneric("continue",function(object,...)standardGeneric("continue"))
setGeneric("pred.mean",function(object,...)standardGeneric("pred.mean"))
setGeneric("pred.var",function(object,...)standardGeneric("pred.var"))
setGeneric("filter.mean",function(object,...)standardGeneric("filter.mean"))
