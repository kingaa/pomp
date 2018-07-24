## 'hitch' takes (as '...') the workhorse specifications (as R functions or
## C snippets, processes these, and hitches them to 'pomp.fun' objects
## suitable for use in the appropriate slots in 'pomp' objects.

hitch <- function (..., templates,
  obsnames, statenames, paramnames, covarnames,
  PACKAGE, globals, cfile, cdir, shlib.args,
  verbose = getOption("verbose", FALSE)) {

  ep <- character(0)

  if (missing(templates))
    stop(ep,sQuote("templates")," must be supplied.",call.=FALSE)

  if (missing(cfile)) cfile <- NULL
  if (missing(cdir)) cdir <- NULL
  if (missing(PACKAGE)) PACKAGE <- NULL
  if (missing(globals)) globals <- NULL
  if (missing(shlib.args)) shlib.args <- NULL
  PACKAGE <- as.character(PACKAGE)

  ## defaults for names of states, parameters, observations, and covariates
  if (missing(statenames)) statenames <- NULL
  if (missing(paramnames)) paramnames <- NULL
  if (missing(obsnames)) obsnames <- NULL
  if (missing(covarnames)) covarnames <- NULL

  statenames <- as.character(statenames)
  paramnames <- as.character(paramnames)
  obsnames <- as.character(obsnames)
  covarnames <- as.character(covarnames)

  if (anyDuplicated(statenames)) {
    stop(ep,"all ",sQuote("statenames")," must be unique", call.=FALSE)
  }
  if (anyDuplicated(paramnames)) {
    stop(ep,"all ",sQuote("paramnames")," must be unique", call.=FALSE)
  }
  if (anyDuplicated(obsnames)) {
    stop(ep,"all ",sQuote("obsnames")," must be unique", call.=FALSE)
  }
  if (anyDuplicated(covarnames)) {
    stop(ep,"all ",sQuote("covarnames")," must be unique", call.=FALSE)
  }

  horses <- list(...)
  snippets <- horses[vapply(horses,is,logical(1),"Csnippet")]
  snippets <- lapply(snippets,slot,"text")

  if (length(snippets) > 0) {
    lib <- tryCatch(
      do.call(
        Cbuilder,
        c(
          list(
            templates=templates,
            dir=cdir,
            name=cfile,
            obsnames=obsnames,
            statenames=statenames,
            paramnames=paramnames,
            covarnames=covarnames,
            globals=globals,
            shlib.args=shlib.args,
            verbose=verbose
          ),
          snippets
        )
      ),
      error = function (e) {
        stop("error in building shared-object library from C snippets: ",
          conditionMessage(e),call.=FALSE)
      }
    )
    libname <- lib$name
  } else {
    libname <- ""
    lib <- NULL
  }

  funs <- vector(mode="list",length=length(horses))
  names(funs) <- names(horses)

  for (s in names(funs)) {
    funs[[s]] <- pomp.fun(
      f=horses[[s]],
      slotname=templates[[s]]$slotname,
      PACKAGE=PACKAGE,
      proto=templates[[s]]$proto,
      Cname=templates[[s]]$Cname,
      libname=libname,
      statenames=statenames,
      paramnames=paramnames,
      obsnames=obsnames,
      covarnames=covarnames
    )
  }

  list(funs=funs,lib=lib)
}

Cbuilder <- function (..., templates, name = NULL, dir = NULL,
  statenames, paramnames, covarnames, obsnames,
  globals, shlib.args = NULL,
  verbose = getOption("verbose",FALSE))
{

  name <- cleanForC(name)
  statenames <- cleanForC(statenames)
  paramnames <- cleanForC(paramnames)
  covarnames <- cleanForC(covarnames)
  obsnames <- cleanForC(obsnames)

  globals <- as(globals,"character")

  snippets <- list(...)

  ## 'registry' holds a list of functions to register
  registry <- c("__pomp_load_stack_incr","__pomp_load_stack_decr")
  ## which utilities are needed?
  utils <- which(sapply(seq_along(pomp_templates$utilities),
    function(x) any(grepl(pomp_templates$utilities[[x]]$trigger,
      snippets))))

  ## rely on "-I" flags under *nix
  if (.Platform$OS.type=="unix") {
    pompheader <- "pomp.h"
  } else {
    pompheader <- system.file("include/pomp.h",package="pomp") # nocov
  }

  ## some information to help make file (and filename) unique
  timestamp <- format(Sys.time(),"%Y-%m-%d %H:%M:%OS3 %z")
  salt <- paste(format(as.hexmode(ceiling(runif(n=4L,max=2^24))),
    upper.case=TRUE),collapse="")

  ## string 'csrc' will hold the full C source code
  csrc <- ""
  out <- textConnection(object="csrc",open="w",local=TRUE)

  cat(file=out,render(pomp_templates$file$header,
    pompheader=pompheader,timestamp=timestamp,salt=salt))

  cat(file=out,globals,"\n\n")

  for (u in utils)
    cat(file=out,pomp_templates$utilities[[u]]$header)

  ## now we write the snippets, using the templates provided
  for (snip in names(snippets)) {
    ## we add each 'Cname' to the 'registry'
    registry <- c(registry,templates[[snip]]$Cname)
    ## define variables
    cat(file=out,render("\n/* C snippet: '{%snip%}' */\n",snip=snip))
    for (k in seq_along(templates[[snip]]$vars)) {
      vtpl <- templates[[snip]]$vars[[k]]
      nm <- eval(vtpl$names)
      for (v in seq_along(nm)) {
        cat(file=out,render(
          pomp_templates$define,
          variable=nm[v],
          cref=render(vtpl$cref,v=v-1)
        ))
      }
    }
    ## render the C snippet in context
    cat(file=out,render(templates[[snip]]$header),
      snippets[[snip]],templates[[snip]]$footer)
    ## undefine variables
    for (k in seq_along(templates[[snip]]$vars)) {
      vtpl <- templates[[snip]]$vars[[k]]
      nm <- eval(vtpl$names)
      for (v in nm) {
        cat(file=out,render(pomp_templates$undefine,variable=v))
      }
    }
  }

  ## load/unload stack handling codes
  cat(file=out,pomp_templates$stackhandling)

  ## registration
  cat(file=out,render(pomp_templates$registration$header))
  for (v in utils)
    cat(file=out,pomp_templates$utilities[[v]]$reg)
  for (v in registry)
    cat(file=out,render(pomp_templates$registration$main,fun=v))
  cat(file=out,pomp_templates$registration$footer)

  close(out)

  csrc <- paste(csrc,collapse="\n")

  if (length(name)==0)
    name <- paste0("pomp_",digest(csrc,serialize=FALSE))

  csrc <- render(csrc,name=name)

  tryCatch(
    pompCompile(
      fname=name,
      direc=srcDir(dir,verbose=verbose),
      src=csrc,
      shlib.args=shlib.args,
      verbose=verbose
    ),
    error = function (e) {
      stop("in ",sQuote("Cbuilder"),": compilation error: ",conditionMessage(e),call.=FALSE)
    }
  )

  invisible(list(name=name,dir=dir,src=csrc))
}

srcDir <- function (dir, verbose) {
  if (is.null(dir)) {
    pid <- Sys.getpid()
    dir <- file.path(tempdir(),pid)
  }
  if (!dir.exists(dir)) {
    if (verbose) cat("creating C snippet directory ",sQuote(dir),"\n")
    tryCatch(
      {
        dir.create(dir,recursive=TRUE,showWarnings=FALSE,mode="0700")
        stopifnot(dir.exists(dir))
      },
      error = function (e) {
        stop("cannot create cache directory ",sQuote(dir),call.=FALSE)
      }
    )
  }
  dir
}

pompCompile <- function (fname, direc, src, shlib.args = NULL, verbose) {

  stem <- file.path(direc,fname)
  if (.Platform$OS.type=="windows")
    stem <- gsub("\\","/",stem,fixed=TRUE) # nocov

  modelfile <- paste0(stem,".c")

  tryCatch(
    cat(src,file=modelfile),
    error = function (e) {
      stop("cannot write file ",sQuote(modelfile),call.=FALSE)   #nocov
    }
  )
  if (verbose) cat("model codes written to",sQuote(modelfile),"\n")

  cflags <- Sys.getenv("PKG_CPPFLAGS")
  cflags <- paste0("PKG_CPPFLAGS=\"",
    if (nchar(cflags)>0) paste0(cflags," ") else "",
    "-I",system.file("include",package="pomp"),
    " -I",getwd(),"\"")

  shlib.args <- as.character(shlib.args)

  solib <- paste0(stem,.Platform$dynlib.ext)
  if (verbose) cat("compiling",sQuote(solib),"\n")
  tryCatch(
    {
      rv <- system2(
        command=R.home("bin/R"),
        args=c("CMD","SHLIB","-c","-o",solib,modelfile,shlib.args),
        env=cflags,
        wait=TRUE,
        stdout=TRUE,
        stderr=TRUE
      )
    },
    error = function (e) {
      stop("error compiling C snippets: ",conditionMessage(e),call.=FALSE) #nocov
    }
  )
  stat <- as.integer(attr(rv,"status"))
  if (length(stat) > 0 && stat != 0L) {
    stop("cannot compile shared-object library ",sQuote(solib),": status = ",stat,
      "\ncompiler messages:\n",paste(rv,collapse="\n"),call.=FALSE)
  } else if (verbose) {
    cat("compiler messages:",rv,sep="\n")
  }

  invisible(solib)
}

cleanForC <- function (text) {
  text <- as.character(text)
  text <- gsub("\\.","_",text)
  text <- gsub("-","_",text)
  text
}

render <- function (template, ...) {
  vars <- list(...)
  if (length(vars)==0) return(template)
  n <- sapply(vars,length)
  if (!all((n==max(n))|(n==1)))
    stop("in ",sQuote("render")," incommensurate lengths of replacements",call.=FALSE)
  short <- which(n==1)
  n <- max(n)
  for (i in short) vars[[i]] <- rep(vars[[i]],n)

  retval <- vector(mode="list",length=n)
  for (i in seq_len(n)) {
    tpl <- template
    for (v in names(vars)) {
      src <- sprintf("\\{%%%s%%\\}",v)
      tgt <- vars[[v]][i]
      tpl <- gsub(src,tgt,tpl)
    }
    retval[[i]] <- tpl
  }
  do.call(paste0,retval)
}

## TEMPLATES

pomp_templates <- list(
  define="#define {%variable%}\t\t({%cref%})\n",
  undefine="#undef {%variable%}\n",
  file=list(
    header="/* pomp C snippet file: {%name%} */\n/* Time: {%timestamp%} */\n/* Salt: {%salt%} */\n\n#include <{%pompheader%}>\n#include <R_ext/Rdynload.h>\n\n"
  ),
  utilities=list(
    periodic_bspline_basis=list(
      trigger="periodic_bspline_basis_eval",
      header="static void (*__pomp_periodic_bspline_basis_eval)(double, double, int, int, double*);\n#define periodic_bspline_basis_eval(X,Y,M,N,Z)\t(__pomp_periodic_bspline_basis_eval((X),(Y),(M),(N),(Z)))
\n",
      reg="__pomp_periodic_bspline_basis_eval = (void (*)(double,double,int,int,double*)) R_GetCCallable(\"pomp\",\"periodic_bspline_basis_eval\");\n"
    ),
    get_pomp_userdata_int=list(
      trigger="get_pomp_userdata_int",
      header="static const int * (*__pomp_get_pomp_userdata_int)(const char *);\n#define get_pomp_userdata_int(X)\t(__pomp_get_pomp_userdata_int(X))\n",
      reg="__pomp_get_pomp_userdata_int = (const int *(*)(const char*)) R_GetCCallable(\"pomp\",\"get_pomp_userdata_int\");\n"
    ),
    get_pomp_userdata_double=list(
      trigger="get_pomp_userdata_double",
      header="static const double * (*__pomp_get_pomp_userdata_double)(const char *);\n#define get_pomp_userdata_double(X)\t(__pomp_get_pomp_userdata_double(X))\n",
      reg="__pomp_get_pomp_userdata_double = (const double *(*)(const char*)) R_GetCCallable(\"pomp\",\"get_pomp_userdata_double\");\n"
    ),
    get_pomp_userdata=list(
      trigger="get_pomp_userdata(\\b|[^_])",
      header="static const SEXP (*__pomp_get_pomp_userdata)(const char *);\n#define get_pomp_userdata(X)\t(__pomp_get_pomp_userdata(X))\n",
      reg="__pomp_get_pomp_userdata = (const SEXP (*)(const char*)) R_GetCCallable(\"pomp\",\"get_pomp_userdata\");\n"
    )
  ),
  stackhandling="\nstatic int __pomp_load_stack = 0;\n\nvoid __pomp_load_stack_incr (void) {++__pomp_load_stack;}\n\nvoid __pomp_load_stack_decr (int *val) {*val = --__pomp_load_stack;}\n",
  registration=list(
    header="\nvoid R_init_{%name%} (DllInfo *info)\n{\n",
    main="R_RegisterCCallable(\"{%name%}\", \"{%fun%}\", (DL_FUNC) {%fun%});\n",
    footer="}\n\n"
  )
)
