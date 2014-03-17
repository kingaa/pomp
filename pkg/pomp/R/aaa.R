## .onAttach <- function (...) {
##   version <- library(help=pomp)$info[[1L]]
##   version <- strsplit(version[pmatch("Version",version)]," ")[[1L]]
##   version <- version[nchar(version)>0][2L]
##   packageStartupMessage("This is pomp version ",version,"\n")
## }

if (!exists("paste0",where="package:base")) {
  paste0 <- function(...) paste(...,sep="")
}
