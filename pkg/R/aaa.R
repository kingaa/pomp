##.onAttach <- function (...) {
##   version <- library(help=pomp)$info[[1]]
##   version <- strsplit(version[pmatch("Version",version)]," ")[[1]]
##   version <- version[nchar(version)>0][2]
##   cat("This is pomp version ",version,"\n\n",sep="")
##   cat("See the NEWS file for important information\n")
##}

setGeneric("print",function(x,...)standardGeneric("print"))
setGeneric("plot",function(x,y,...)standardGeneric("plot"))
setGeneric("summary",function(object,...)standardGeneric("summary"))
setGeneric("simulate",function(object,nsim=1,seed=NULL,...)standardGeneric("simulate"))
setGeneric("time",function(x,...)standardGeneric("time"))
setGeneric("coef",function(object,...)standardGeneric("coef"))
setGeneric("logLik",function(object,...)standardGeneric("logLik"))
setGeneric("window",function(x,...)standardGeneric("window"))
setGeneric("continue",function(object,...)standardGeneric("continue"))
