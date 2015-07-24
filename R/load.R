pompLoad.internal <- function (object, ..., verbose = getOption("verbose",FALSE)) {
  for (lib in object@solibs) {
    if (!is.loaded("__pomp_load_stack_incr",PACKAGE=lib$name)) {
      dir <- pompSrcDir(lib$dir)
      solib <- file.path(dir,paste0(lib$name,.Platform$dynlib.ext))
      if (file.exists(solib)) {
        dyn.load(solib)
      } else {
        pompCompile(fname=lib$name,direc=dir,src=lib$src,verbose=verbose)
        dyn.load(solib)
      }
      if (verbose) cat("loading",sQuote(solib),"\n")
    }
    .Call(load_stack_incr,lib$name)
  }
  invisible(NULL)
}
 
pompUnload.internal <- function (object, ..., verbose = getOption("verbose",FALSE)) {
  for (lib in object@solibs) {
    if (is.loaded("__pomp_load_stack_decr",PACKAGE=lib$name)) {
      st <- .Call(load_stack_decr,lib$name)
      if (st==0) {
        dir <- pompSrcDir(lib$dir)
        solib <- file.path(dir,paste0(lib$name,.Platform$dynlib.ext))
        dyn.unload(solib)
        if (verbose) cat("unloading",sQuote(solib),"\n")
      }
    }
  }
  invisible(NULL)
}

setMethod("pompLoad",
          signature=signature(object='pomp'),
          definition = function (object, ...) {
            pompLoad.internal(object,...)
          })

setMethod("pompUnload",
          signature=signature(object='pomp'),
          definition = function (object, ...) {
            pompUnload.internal(object,...)
          })
