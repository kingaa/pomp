.onAttach <- function (...) {
  version <- library(help=pomp)$info[[1]]
  version <- strsplit(version[pmatch("Version",version)]," ")[[1]]
  version <- version[nchar(version)>0][2]
  cat("This is pomp version ",version,"\n\n",sep="")
  paste(
        "IMPORTANT NOTICE:\n",
        "The default behaviors of ",sQuote("simulate")," and ",sQuote("trajectory"),
        " have changed as of release 0.34-1. ",
        "You can ensure that your code will continue to function as you intend by specifying the values of the ",
        sQuote("times")," and ",sQuote("t0"),
        " arguments to these functions, thus removing dependence of your code on the defaults. ",
        "In the meantime, using ",sQuote("simulate")," or ",sQuote("trajectory"),
        " in such a way as to rely on the default will produce a warning. ",
        "These warnings will be removed in a future release.\n\n",
        "See the documentation (",dQuote("pomp?simulate"),", ",dQuote("pomp?trajectory"),
        ") for a description of the new default behaviors.\n\n",
        "Subscribe to the pomp-announce list (go to pomp.r-forge.r-project.org) ",
        "to receive email notification about new releases of pomp including ",
        "descriptions of feature additions, changes, and fixes.\n\n",
        sep=""
        ) -> msg
  cat(msg,fill=TRUE)
}
