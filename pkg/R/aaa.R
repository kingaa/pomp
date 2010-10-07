##.onAttach <- function (...) {
##   version <- library(help=pomp)$info[[1]]
##   version <- strsplit(version[pmatch("Version",version)]," ")[[1]]
##   version <- version[nchar(version)>0][2]
##   cat("This is pomp version ",version,"\n\n",sep="")
##   cat("See the NEWS file for important information\n")
##}
