pompLoad.internal <- function (object, ..., verbose = getOption("verbose",FALSE)) {
  for (lib in object@solibfile) {
    if (!is.loaded("__pomp_load_stack_incr",PACKAGE=lib[1])) {
      dyn.load(lib[2])
      if (verbose) cat("loading",sQuote(lib[2]),"\n")
    }
    .Call(load_stack_incr,lib[1])
  }
  invisible(NULL)
}
 
pompUnload.internal <- function (object, ..., verbose = getOption("verbose",FALSE)) {
  for (lib in object@solibfile) {
    if (is.loaded("__pomp_load_stack_decr",PACKAGE=lib[1])) {
      st <- .Call(load_stack_decr,lib[1])
      if (st==0) {
        dyn.unload(lib[2])
        if (verbose) cat("unloading",sQuote(lib[2]),"\n")
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
