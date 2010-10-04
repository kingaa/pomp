.onAttach <- function(...) {
  version <- packageDescription(pkg="pomp",fields="Version")
  cat("This is pomp version ",version,"\n",sep="")
  cat("\n",sep="")
}
