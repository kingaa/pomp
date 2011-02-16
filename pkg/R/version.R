version <- function (at.least = NULL) {
  version <- library(help=pomp)$info[[1]]
  version <- strsplit(version[pmatch("Version",version)]," ")[[1]]
  version <- version[nchar(version)>0][2]
  version <- as.numeric(strsplit(version,"[-.]")[[1]])
  if (is.null(at.least)) {
    list(major=version[1],minor=version[2],rev=version[3])
  } else {
    minv <- as.numeric(strsplit(as.character(at.least),"[-.]")[[1]])
    (version[1]>minv[1]) ||
    (version[1]==minv[1]) && (version[2]>minv[2]) ||
    (version[1]==minv[1]) && (version[2]==minv[2]) && (version[3]>=minv[3])
  }
}

