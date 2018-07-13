pompCBuilder <- function (name = NULL, dir = NULL,
  statenames, paramnames, covarnames, obsnames,
  ..., globals, shlib.args = NULL,
  templates = snippet_templates,
  verbose = getOption("verbose",FALSE))
{

  if (!is.null(name)) name <- cleanForC(name)

  if (is(globals,"Csnippet")) globals <- globals@text

  statenames <- cleanForC(statenames)
  paramnames <- cleanForC(paramnames)
  covarnames <- cleanForC(covarnames)
  obsnames <- cleanForC(obsnames)

  snippets <- list(...)

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

  ## list of functions to register
  registry <- c("__pomp_load_stack_incr","__pomp_load_stack_decr")

  ## string 'csrc' will hold the full C source code
  csrc <- ""
  out <- textConnection(object="csrc",open="w",local=TRUE)

  cat(file=out,render(pomp_templates$file$header,
    pompheader=pompheader,timestamp=timestamp,salt=salt))

  cat(file=out,globals,"\n\n")

  needsmap <- list(
    statenames=statenames,
    derivs=statenames,
    before_n_after=statenames,
    paramnames=paramnames,
    transforms=paramnames,
    covarnames=covarnames,
    obsnames=obsnames,
    lik=TRUE,
    loglik=TRUE
  )

  ## now we write the snippets, using the templates provided
  ## we add each 'Cname' to the 'registry'
  for (snip in names(snippets)) {
    registry <- c(registry,templates[[snip]]$Cname)
    do.call(defmacros,c(list(out),needsmap[templates[[snip]]$needs]))
    cat(file=out,render(templates[[snip]]$header),
      callable.decl(snippets[[snip]]),snippets[[snip]],
      templates[[snip]]$footer)
    do.call(undefmacros,c(list(out),needsmap[templates[[snip]]$needs]))
  }

  ## load/unload stack handling codes
  cat(file=out,pomp_templates$stackhandling)

  ## registration
  cat(file=out,render(pomp_templates$registration$header))
  for (v in registry)
    cat(file=out,render(pomp_templates$registration$main,fun=v))
  cat(file=out,pomp_templates$registration$footer)

  close(out)

  csrc <- paste(csrc,collapse="\n")

  if (is.null(name))
    name <- paste0("pomp_",digest(csrc,serialize=FALSE))

  csrc <- render(csrc,name=name)

  tryCatch(
    pompCompile(
      fname=name,
      direc=pompSrcDir(dir,verbose=verbose),
      src=csrc,
      shlib.args=shlib.args,
      verbose=verbose
    ),
    error = function (e) {
      stop("in ",sQuote("pompCBuilder"),": compilation error: ",conditionMessage(e),call.=FALSE)
    }
  )

  invisible(list(name=name,dir=dir,src=csrc))
}

defmacros <- function (out, covarnames, paramnames, statenames, obsnames,
  before_n_after, derivs, transforms, lik = FALSE, loglik = FALSE) {
  f1 <- function (nm, ptr, idx) {
    for (v in seq_along(nm)) {
      cat(file=out,render(pomp_templates$define$var,variable=nm[v],ptr=ptr,ilist=idx,
        index=as.integer(v-1)))
    }
  }

  if (!missing(covarnames)) f1(covarnames,"__covars","__covindex")
  if (!missing(paramnames)) f1(paramnames,"__p","__parindex")
  if (!missing(statenames)) f1(statenames,"__x","__stateindex")
  if (!missing(obsnames)) f1(obsnames,"__y","__obsindex")
  if (!missing(before_n_after)) {
    f1(paste0(before_n_after,"_1"),"__x1","__stateindex")
    f1(paste0(before_n_after,"_2"),"__x2","__stateindex")
  }
  if (!missing(derivs)) f1(paste0("D",derivs),"__f","__stateindex")
  if (!missing(transforms)) f1(paste0("T",transforms),"__pt","__parindex")

  f2 <- function (nm, ptr) {
    cat(file=out,render(pomp_templates$define$var.alt,variable=nm,ptr=ptr,index=0L))
  }

  if (lik) f2("lik","__lik")
  if (loglik) f2("loglik","__loglik")

  invisible(NULL)
}

undefmacros <- function (out, covarnames, paramnames, statenames, obsnames,
  before_n_after, derivs, transforms, lik = FALSE, loglik = FALSE) {

  f1 <- function (nm) {
    for (v in nm) {
      cat(file=out,render(pomp_templates$undefine$var,variable=v))
    }
  }

  if (!missing(covarnames)) f1(covarnames)
  if (!missing(paramnames)) f1(paramnames)
  if (!missing(statenames)) f1(statenames)
  if (!missing(obsnames)) f1(obsnames)
  if (!missing(before_n_after)) {
    f1(paste0(before_n_after,"_1"))
    f1(paste0(before_n_after,"_2"))
  }
  if (!missing(derivs)) f1(paste0("D",derivs))
  if (!missing(transforms)) f1(paste0("T",transforms))

  if (lik) f1("lik")
  if (loglik) f1("loglik")

  invisible(NULL)
}

pompSrcDir <- function (dir, verbose) {
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
    "-I",system.file("include",package="pomp"),"\"")

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

callable.decl <- function (code) {

  utilities <- list(
    periodic_bspline_basis_eval="\tvoid (*periodic_bspline_basis_eval)(double,double,int,int,double*);\nperiodic_bspline_basis_eval = (void (*)(double,double,int,int,double*)) R_GetCCallable(\"pomp\",\"periodic_bspline_basis_eval\");\n",
    get_pomp_userdata_int="\tconst int * (*get_pomp_userdata_int)(const char *);\nget_pomp_userdata_int = (const int *(*)(const char*)) R_GetCCallable(\"pomp\",\"get_pomp_userdata_int\");\n",
    get_pomp_userdata_double="\tconst double * (*get_pomp_userdata_double)(const char *);\nget_pomp_userdata_double = (const double *(*)(const char*)) R_GetCCallable(\"pomp\",\"get_pomp_userdata_double\");\n",
    `get_pomp_userdata(\\b|[^_])`="\tconst SEXP (*get_pomp_userdata)(const char *);\nget_pomp_userdata = (const SEXP (*)(const char*)) R_GetCCallable(\"pomp\",\"get_pomp_userdata\");\n"
  )

  fns <- vapply(names(utilities),grepl,logical(1),code,perl=TRUE)
  do.call(paste0,utilities[fns])
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
  define=list(
    var="#define {%variable%}\t({%ptr%}[{%ilist%}[{%index%}]])\n",
    var.alt="#define {%variable%}\t({%ptr%}[{%index%}])\n"
  ),
  undefine=list(
    var="#undef {%variable%}\n"
  ),
  file=list(
    header="/* pomp C snippet file: {%name%} */\n/* Time: {%timestamp%} */\n/* Salt: {%salt%} */\n\n#include <{%pompheader%}>\n#include <R_ext/Rdynload.h>\n\n"
  ),
  stackhandling="\nstatic int __pomp_load_stack = 0;\n\nvoid __pomp_load_stack_incr (void) {++__pomp_load_stack;}\n\nvoid __pomp_load_stack_decr (int *val) {*val = --__pomp_load_stack;}\n",
  registration=list(
    header="\nvoid R_init_{%name%} (DllInfo *info)\n{\n",
    main="R_RegisterCCallable(\"{%name%}\", \"{%fun%}\", (DL_FUNC) {%fun%});\n",
    footer="}\n\n"
  )
)
