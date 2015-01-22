pompExample <- function (example, ..., envir = .GlobalEnv) {
  example <- as.character(substitute(example))
  pomp.dir <- system.file("examples",package="pomp")
  exampleDirs <- getOption("pomp.examples",default=pomp.dir)
  names(exampleDirs) <- exampleDirs
  if (example=="") {
    avlbl <- lapply(exampleDirs,list.files,pattern=".+?R$")
    avlbl <- lapply(avlbl,function(x) gsub("\\.R$","",x))
    for (dir in exampleDirs) {
      cat("examples in ",dir,":\n",sep="")
      print(avlbl[[dir]])
    }
  } else {
## the following needed from R/3.0.0 to R/3.1.1:
    dots <- list(...)
    evalEnv <- if (length(dots)>0) list2env(list(...)) else new.env()
## the following will work from R/3.1.2
##  evalEnv <- list2env(list(...))
    file <- c(lapply(exampleDirs,list.files,
                     pattern=paste0(example,".R"),
                     full.names=TRUE),
              recursive=TRUE)
    if (length(file)>1) {
      warning("using ",sQuote(file[1])," from ",sQuote(names(file)[1]))
    }
    objs <- source(file[1],local=evalEnv)
    if (is.null(envir)) {
      obj <- setNames(lapply(objs$value,get,envir=evalEnv),objs$value)
    } else if (is.environment(envir)) {
      for (i in seq_along(objs$value)) {
        assign(objs$value[i],
               get(objs$value[i],envir=evalEnv),
               envir=envir)
      }
      cat("newly created object(s):\n",objs$value,"\n")
      obj <- NULL
    } else {
      stop(sQuote("envir")," must be an environment or NULL")
    }
    invisible(obj)
  }
}
