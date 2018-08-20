##' Examples of the construction of POMP models
##'
##' \code{pompExample} loads pre-built example \sQuote{pomp} objects.
##'
##' Directories listed in the global option \code{pomp.examples} (which can be
##' changed using \code{options()}) are searched for file named
##' \file{<example>.R}.  If found, this file will be \code{source}d in a
##' temporary environment.  Additional arguments to \code{pompExample} define
##' variables within this environment and will therefore be available when the
##' code in \file{<example>.R} is \code{source}d.
##'
##' The codes that construct these \sQuote{pomp} objects can be found in the
##' \file{examples} directory in the installed package.  Do
##' \code{system.file("examples",package="pomp"))} to find this directory.
##'
##' @name pompExample
##' @rdname pompExample
##' @include pomp.R simulate.R
##' @importFrom stats setNames
##' @author Aaron A. King
##' @keywords models datasets
##' @family pomp examples
##'
##' @param example example to load given as a name or literal character string.
##' Evoked without an argument, \code{pompExample} lists all available
##' examples.
##' @param \dots additional arguments define symbols in the environment within
##' which the example code is executed.
##' @param show logical; if \code{TRUE}, display, but do not execute, the
##' example code.
##' @param envir the environment into which the objects should be loaded.  If
##' \code{envir=NULL}, then the created objects are returned in a list.
##'
##' @return By default, \code{pompExample} has the side effect of creating one
##' or more objects in the global workspace.  If \code{envir=NULL}, there are
##' no side effects; rather, the objects are returned as a list.
##'
##' @examples
##'
##'   pompExample()
##'   pompExample(sir)
##'   pompExample("gompertz")
##'   pompExample(ricker,envir=NULL)
##' \dontrun{
##'   pompExample(bbs,show=TRUE)
##' }
##'
NULL

##' @export
pompExample <- function (example, ..., show = FALSE, envir = .GlobalEnv) {
  example <- as.character(substitute(example))
  ep <- paste0("in ",sQuote("pompExample"),": ")
  pomp.dir <- system.file("examples",package="pomp")
  exampleDirs <- getOption("pomp.examples",default=pomp.dir)
  names(exampleDirs) <- exampleDirs
  show <- as.logical(show)
  if (example=="") {
    avlbl <- lapply(exampleDirs,list.files,pattern=".+?R$")
    avlbl <- lapply(avlbl,function(x) gsub("\\.R$","",x))
    for (dir in exampleDirs) {
      cat("examples in ",dir,":\n",sep="")
      print(avlbl[[dir]])
    }
  } else {
    evalEnv <- list2env(list(...))
    file <- c(lapply(exampleDirs,list.files,
      pattern=paste0(example,".R"),
      full.names=TRUE),
      recursive=TRUE)
    if (length(file)<1) {
      stop(ep,"cannot find file ",
        sQuote(paste0(example,".R")),call.=FALSE)
    }
    if (length(file)>1) {
      warning(ep,"using ",sQuote(file[1])," from ",sQuote(names(file)[1]),call.=FALSE)
    }
    if (show) {
      file.show(file[1])
      return(invisible(NULL))
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
      stop(ep,sQuote("envir")," must be an environment or NULL",call.=FALSE)
    }
    invisible(obj)
  }
}
