pompCBuilder <- function (name = NULL, dir = NULL,
  statenames, paramnames, covarnames, obsnames,
  ..., globals, shlib.args = NULL,
  templates = pomp_templates,
  verbose = getOption("verbose",FALSE))
{

  if (!is.null(name)) name <- cleanForC(name)

  if (is(globals,"Csnippet")) globals <- globals@text

  statenames <- cleanForC(statenames)
  paramnames <- cleanForC(paramnames)
  covarnames <- cleanForC(covarnames)
  obsnames <- cleanForC(obsnames)

  snippets <- list(...)
  snipnms <- names(snippets)

  has <- function (...) any(c(...) %in% snipnms)

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

  cat(file=out,render(templates$file$header,
    pompheader=pompheader,timestamp=timestamp,salt=salt))

  ##    for (f in utility.fns) {
  ##        cat(file=out,f)
  ##    }

  cat(file=out,globals,"\n\n")

  ## state-variable/parameter/covariates macros
  for (v in seq_along(paramnames)) {
    cat(file=out,render(templates$define$var,variable=paramnames[v],ptr='__p',ilist='__parindex',index=as.integer(v-1)))
  }
  for (v in seq_along(covarnames)) {
    cat(file=out,render(templates$define$var,variable=covarnames[v],ptr='__covars',ilist='__covindex',index=as.integer(v-1)))
  }

  for (v in seq_along(statenames)) {
    cat(file=out,render(templates$define$var,variable=statenames[v],ptr='__x',ilist='__stateindex',index=as.integer(v-1)))
  }

  ## initializer function
  if (has("initializer")) {
    registry <- c(registry,fnames[["initializer"]])
    cat(file=out,render(templates$initializer$header),
      callable.decl(snippets$initializer),snippets$initializer,
      templates$initializer$footer)
  }

  ## Euler step function
  if (has("step.fn")) {
    registry <- c(registry,fnames[["step.fn"]])
    cat(file=out,render(templates$step.fn$header),
      callable.decl(snippets$step.fn),snippets$step.fn,
      templates$step.fn$footer)
  }

  ## Gillespie rate function
  if (has("rate.fn")) {
    registry <- c(registry,fnames[["rate.fn"]])
    cat(file=out,render(templates$rate.fn$header),
      callable.decl(snippets$rate.fn),snippets$rate.fn,
      templates$rate.fn$footer)
  }

  if (has("rmeasure","dmeasure")) {

    for (v in seq_along(obsnames)) {
      cat(file=out,render(templates$define$var,variable=obsnames[v],ptr='__y',ilist='__obsindex',index=as.integer(v-1)))
    }

    ## rmeasure function
    if (has("rmeasure")) {
      registry <- c(registry,fnames[["rmeasure"]])
      cat(file=out,render(templates$rmeasure$header),
        callable.decl(snippets$rmeasure),snippets$rmeasure,
        templates$rmeasure$footer)
    }

    ## dmeasure function
    if (has("dmeasure")) {
      cat(file=out,render(templates$define$var.alt,variable="lik",ptr='__lik',index=0L))
      registry <- c(registry,fnames[["dmeasure"]])
      cat(file=out,render(templates$dmeasure$header),
        callable.decl(snippets$dmeasure),snippets$dmeasure,
        templates$dmeasure$footer)
      cat(file=out,render(templates$undefine$var,variable="lik"))
    }

    for (v in obsnames) {
      cat(file=out,render(templates$undefine$var,variable=v))
    }
  }

  ## rprior function
  if (has("rprior")) {
    registry <- c(registry,fnames[["rprior"]])
    cat(file=out,render(templates$rprior$header),
      callable.decl(snippets$rprior),snippets$rprior,
      templates$rprior$footer)
  }

  ## dprior function
  if (has("dprior")) {
    cat(file=out,render(templates$define$var.alt,variable="lik",ptr='__lik',index=0L))
    registry <- c(registry,fnames[["dprior"]])
    cat(file=out,render(templates$dprior$header),
      callable.decl(snippets$dprior),snippets$dprior,
      templates$dprior$footer)
    cat(file=out,render(templates$undefine$var,variable="lik"))
  }

  if (has("fromEstimationScale","toEstimationScale")) {

    for (v in seq_along(paramnames)) {
      cat(file=out,render(templates$define$var,variable=paste0("T",paramnames[v]),
        ptr='__pt',ilist='__parindex',index=as.integer(v-1)))
    }

    ## parameter transformation function
    if (has("fromEstimationScale")) {
      registry <- c(registry,fnames[["fromEstimationScale"]])
      cat(file=out,render(templates$fromEstimationScale$header),
        callable.decl(snippets$fromEstimationScale),snippets$fromEstimationScale,
        templates$toEstimationScale$footer)
    }

    ## inverse parameter transformation function
    if (has("toEstimationScale")) {
      registry <- c(registry,fnames[["toEstimationScale"]])
      cat(file=out,render(templates$toEstimationScale$header),
        callable.decl(snippets$toEstimationScale),snippets$toEstimationScale,
        templates$toEstimationScale$footer)
    }

    for (v in paramnames) {
      cat(file=out,render(templates$undefine$var,variable=paste0("T",v)))
    }
  }

  ## skeleton function
  if (has("skeleton")) {
    registry <- c(registry,fnames[["skeleton"]])
    for (v in seq_along(statenames)) {
      cat(file=out,render(templates$define$var,variable=paste0("D",statenames[v]),
        ptr='__f',ilist='__stateindex',index=as.integer(v-1)))
    }
    cat(file=out,render(templates$skeleton$header),
      callable.decl(snippets$skeleton),snippets$skeleton,
      templates$skeleton$footer)
    for (v in statenames) {
      cat(file=out,render(templates$undefine$var,variable=paste0("D",v)))
    }
  }

  ## undefine state variables
  for (v in statenames) {
    cat(file=out,render(templates$undefine$var,variable=v))
  }

  if (has("dens.fn")) {
    cat(file=out,render(templates$define$var.alt,variable="loglik",ptr='__loglik',index=0L))
    for (v in seq_along(statenames)) {
      cat(file=out,render(templates$define$var,variable=paste0(statenames[v],"_1"),ptr='__x1',ilist='__stateindex',index=as.integer(v-1)))
      cat(file=out,render(templates$define$var,variable=paste0(statenames[v],"_2"),ptr='__x2',ilist='__stateindex',index=as.integer(v-1)))
    }

    registry <- c(registry,fnames[["dens.fn"]])
    cat(file=out,render(templates$dens.fn$header),
      callable.decl(snippets$dens.fn),snippets$dens.fn,
      templates$dens.fn$footer)

    for (v in statenames) {
      cat(file=out,render(templates$undefine$var,variable=paste0(v,"_1")))
      cat(file=out,render(templates$undefine$var,variable=paste0(v,"_2")))
    }
    cat(file=out,render(templates$undefine$var,variable="loglik"))
  }

  ## more undefines
  for (v in paramnames) {
    cat(file=out,render(templates$undefine$var,variable=v))
  }
  for (v in covarnames) {
    cat(file=out,render(templates$undefine$var,variable=v))
  }

  ## load/unload stack handling codes
  cat(file=out,templates$stackhandling)

  ## registration
  cat(file=out,render(templates$registration$header))
  for (v in registry)
    cat(file=out,render(templates$registration$main,fun=v))
  cat(file=out,templates$registration$footer)

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

## the 'fnames' are needed in 'pomp_fun.R'
fnames <- list(
  rmeasure="__pomp_rmeasure",
  dmeasure= "__pomp_dmeasure",
  step.fn="__pomp_stepfn",
  rate.fn="__pomp_ratefn",
  dens.fn="__pomp_densfn",
  skeleton="__pomp_skelfn",
  initializer="__pomp_initializer",
  fromEstimationScale="__pomp_par_trans",
  toEstimationScale="__pomp_par_untrans",
  rprior="__pomp_rprior",
  dprior="__pomp_dprior"
)

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
  stackhandling="
static int __pomp_load_stack = 0;

void __pomp_load_stack_incr (void) {++__pomp_load_stack;}

void __pomp_load_stack_decr (int *val) {*val = --__pomp_load_stack;}
",
  initializer=list(
    header="\nvoid __pomp_initializer (double *__x, const double *__p, double t, const int *__stateindex, const int *__parindex, const int *__covindex, const double *__covars)\n{\n",
    footer="\n}\n\n"
  ),
  rmeasure=list(
    header="\nvoid __pomp_rmeasure (double *__y, const double *__x, const double *__p, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, int __ncovars, const double *__covars, double t)\n{\n",
    footer="\n}\n\n"
  ),
  dmeasure=list(
    header="\nvoid __pomp_dmeasure (double *__lik, const double *__y, const double *__x, const double *__p, int give_log, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, int __ncovars, const double *__covars, double t)\n{\n",
    footer="\n}\n\n"
  ),
  step.fn=list(
    header="\nvoid __pomp_stepfn (double *__x, const double *__p, const int *__stateindex, const int *__parindex, const int *__covindex, int __covdim, const double *__covars, double t, double dt)\n{\n",
    footer="\n}\n\n"
  ),
  rate.fn=list(
    header="\ndouble __pomp_ratefn (int j, double t, double *__x, const double *__p, const int *__stateindex, const int *__parindex, const int *__covindex, int __covdim, const double *__covars)\n{\n  double rate = 0.0;  \n",
    footer="  return rate;\n}\n\n"
  ),
  dens.fn=list(
    header="\nvoid __pomp_densfn (double *__loglik, const double *__x1, const double *__x2, double t1, double t2, const double *__p, const int *__stateindex, const int *__parindex, const int *__covindex, int __ncovars, const double *__covars)\n{\n",
    footer="\n}\n\n"
  ),
  skeleton=list(
    header="\nvoid __pomp_skelfn (double *__f, const double *__x, const double *__p, const int *__stateindex, const int *__parindex, const int *__covindex, int __ncovars, const double *__covars, double t)\n{\n",
    footer="\n}\n\n"
  ),
  fromEstimationScale=list(
    header="\nvoid __pomp_par_trans (double *__pt, const double *__p, const int *__parindex)\n{\n",
    footer="\n}\n\n"
  ),
  toEstimationScale=list(
    header="\nvoid __pomp_par_untrans (double *__pt, const double *__p, const int *__parindex)\n{\n",
    footer="\n}\n\n"
  ),
  rprior=list(
    header="\nvoid __pomp_rprior (double *__p, const int *__parindex)\n{\n",
    footer="\n}\n\n"
  ),
  dprior=list(
    header="\nvoid __pomp_dprior (double *__lik, const double *__p, int give_log, const int *__parindex)\n{\n",
    footer="\n}\n\n"
  ),
  registration=list(
    header="\nvoid R_init_{%name%} (DllInfo *info)\n{\n",
    main="R_RegisterCCallable(\"{%name%}\", \"{%fun%}\", (DL_FUNC) {%fun%});\n",
    footer="}\n\n"
  )
)

## utility.fns <- list()
