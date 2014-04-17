version <- function (at.least = NULL) {
  version <- library(help=mif2)$info[[1]]
  version <- strsplit(version[pmatch("Version",version)]," ")[[1]]
  version <- version[nchar(version)>0][2]
  splv <- as.numeric(strsplit(version,"[-.]")[[1]])
  if (is.null(at.least)) {
    list(major=splv[1],minor=splv[2],rev=splv[3],version.string=version)
  } else {
    minv <- as.numeric(strsplit(as.character(at.least),"[-.]")[[1]])
    (splv[1]>minv[1]) ||
    (splv[1]==minv[1]) && (splv[2]>minv[2]) ||
    (splv[1]==minv[1]) && (splv[2]==minv[2]) && (splv[3]>=minv[3])
  }
}
