pompCBuilder <- function (name = NULL, dir = NULL,
  statenames, paramnames, covarnames, obsnames,
  rmeasure, dmeasure, step.fn, rate.fn, dens.fn, skeleton,
  fromEstimationScale, toEstimationScale,
  initializer, rprior, dprior, globals, shlib.args = NULL,
  verbose = getOption("verbose",FALSE))
{

  if (!is.null(name)) name <- cleanForC(name)

  if (is(globals,"Csnippet")) globals <- globals@text

  statenames <- cleanForC(statenames)
  paramnames <- cleanForC(paramnames)
  covarnames <- cleanForC(covarnames)
  obsnames <- cleanForC(obsnames)

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
  registry <- c("load_stack_incr","load_stack_decr")

  ## string 'csrc' will hold the full C source code
  csrc <- ""
  out <- textConnection(object="csrc",open="w",local=TRUE)

  cat(file=out,render(header$file,pompheader=pompheader,timestamp=timestamp,salt=salt))

  ##    for (f in utility.fns) {
  ##        cat(file=out,f)
  ##    }

  cat(file=out,globals,footer$globals)

  ## state-variable/parameter/covariates macros
  for (v in seq_along(paramnames)) {
    cat(file=out,render(define$var,variable=paramnames[v],ptr='__p',ilist='__parindex',index=as.integer(v-1)))
  }
  for (v in seq_along(covarnames)) {
    cat(file=out,render(define$var,variable=covarnames[v],ptr='__covars',ilist='__covindex',index=as.integer(v-1)))
  }
  for (v in seq_along(statenames)) {
    cat(file=out,render(define$var,variable=statenames[v],ptr='__x',ilist='__stateindex',index=as.integer(v-1)))
  }

  ## initializer function
  if (!missing(initializer)) {
    registry <- c(registry,"initializer")
    cat(file=out,render(header$initializer))
    cat(file=out,callable.decl(initializer))
    cat(file=out,initializer,footer$initializer)
  }

  ## Euler step function
  if (!missing(step.fn)) {
    registry <- c(registry,"stepfn")
    cat(file=out,render(header$step.fn))
    cat(file=out,callable.decl(step.fn))
    cat(file=out,step.fn,footer$step.fn)
  }

  ## Gillespie rate function
  if (!missing(rate.fn)) {
    registry <- c(registry,"ratefn")
    cat(file=out,render(header$rate.fn))
    cat(file=out,callable.decl(rate.fn))
    cat(file=out,rate.fn,footer$rate.fn)
  }

  if (!(missing(rmeasure) && missing(dmeasure))) {

    for (v in seq_along(obsnames)) {
      cat(file=out,render(define$var,variable=obsnames[v],ptr='__y',ilist='__obsindex',index=as.integer(v-1)))
    }

    ## rmeasure function
    if (!missing(rmeasure)) {
      registry <- c(registry,"rmeasure")
      cat(file=out,render(header$rmeasure))
      cat(file=out,callable.decl(rmeasure))
      cat(file=out,rmeasure,footer$rmeasure)
    }

    ## dmeasure function
    if (!missing(dmeasure)) {
      cat(file=out,render(define$var.alt,variable="lik",ptr='__lik',index=0L))
      registry <- c(registry,"dmeasure")
      cat(file=out,render(header$dmeasure))
      cat(file=out,callable.decl(dmeasure))
      cat(file=out,dmeasure,footer$dmeasure)
      cat(file=out,render(undefine$var,variable="lik"))
    }

    for (v in obsnames) {
      cat(file=out,render(undefine$var,variable=v))
    }
  }

  ## rprior function
  if (!missing(rprior)) {
    registry <- c(registry,"rprior")
    cat(file=out,render(header$rprior))
    cat(file=out,callable.decl(rprior))
    cat(file=out,rprior,footer$rprior)
  }

  ## dprior function
  if (!missing(dprior)) {
    cat(file=out,render(define$var.alt,variable="lik",ptr='__lik',index=0L))
    registry <- c(registry,"dprior")
    cat(file=out,render(header$dprior))
    cat(file=out,callable.decl(dprior))
    cat(file=out,dprior,footer$dprior)
    cat(file=out,render(undefine$var,variable="lik"))
  }

  if (!(missing(fromEstimationScale) && missing(toEstimationScale))) {
    for (v in seq_along(paramnames)) {
      cat(file=out,render(define$var,variable=paste0("T",paramnames[v]),
        ptr='__pt',ilist='__parindex',index=as.integer(v-1)))
    }

    ## parameter transformation function
    if (!missing(fromEstimationScale)) {
      registry <- c(registry,"par_trans")
      cat(file=out,render(header$fromEstimationScale))
      cat(file=out,callable.decl(fromEstimationScale))
      cat(file=out,fromEstimationScale,footer$fromEstimationScale)
    }

    ## inverse parameter transformation function
    if (!missing(toEstimationScale)) {
      registry <- c(registry,"par_untrans")
      cat(file=out,render(header$toEstimationScale))
      cat(file=out,callable.decl(toEstimationScale))
      cat(file=out,toEstimationScale,footer$toEstimationScale)
    }

    for (v in paramnames) {
      cat(file=out,render(undefine$var,variable=paste0("T",v)))
    }
  }

  ## skeleton function
  if (!missing(skeleton)) {
    registry <- c(registry,"skelfn")
    for (v in seq_along(statenames)) {
      cat(file=out,render(define$var,variable=paste0("D",statenames[v]),ptr='__f',ilist='__stateindex',index=as.integer(v-1)))
    }
    cat(file=out,render(header$skeleton))
    cat(file=out,callable.decl(skeleton))
    cat(file=out,skeleton,footer$skeleton)
    for (v in statenames) {
      cat(file=out,render(undefine$var,variable=paste0("D",v)))
    }
  }

  ## undefine state variables
  for (v in statenames) {
    cat(file=out,render(undefine$var,variable=v))
  }

  if (!missing(dens.fn)) {
    cat(file=out,render(define$var.alt,variable="loglik",ptr='__loglik',index=0L))
    for (v in seq_along(statenames)) {
      cat(file=out,render(define$var,variable=paste0(statenames[v],"_1"),ptr='__x1',ilist='__stateindex',index=as.integer(v-1)))
      cat(file=out,render(define$var,variable=paste0(statenames[v],"_2"),ptr='__x2',ilist='__stateindex',index=as.integer(v-1)))
    }

    registry <- c(registry,"densfn")
    cat(file=out,render(header$dens.fn))
    cat(file=out,callable.decl(dens.fn))
    cat(file=out,dens.fn,footer$dens.fn)

    for (v in statenames) {
      cat(file=out,render(undefine$var,variable=paste0(v,"_1")))
      cat(file=out,render(undefine$var,variable=paste0(v,"_2")))
    }
    cat(file=out,render(undefine$var,variable="loglik"))
  }

  ## more undefines
  for (v in paramnames) {
    cat(file=out,render(undefine$var,variable=v))
  }
  for (v in covarnames) {
    cat(file=out,render(undefine$var,variable=v))
  }

  ## load/unload stack handling codes
  cat(file=out,stackhandling)

  ## registration
  cat(file=out,render(header$registration))
  for (v in registry) cat(file=out,render(registration,fun=v))
  cat(file=out,footer$registration)

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
  fns <- vapply(names(decl),grepl,logical(1),code,perl=TRUE)
  do.call(paste0,decl[fns])
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

define <- list(
  var="#define {%variable%}\t({%ptr%}[{%ilist%}[{%index%}]])\n",
  var.alt="#define {%variable%}\t({%ptr%}[{%index%}])\n"
)

undefine <- list(
  var="#undef {%variable%}\n"
)

header <- list(
  file="/* pomp C snippet file: {%name%} */\n/* Time: {%timestamp%} */\n/* Salt: {%salt%} */\n\n#include <{%pompheader%}>\n#include <R_ext/Rdynload.h>\n\n",
  registration="\nvoid R_init_{%name%} (DllInfo *info)\n{\n",
  rmeasure="\nvoid __pomp_rmeasure (double *__y, const double *__x, const double *__p, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, int __ncovars, const double *__covars, double t)\n{\n",
  dmeasure= "\nvoid __pomp_dmeasure (double *__lik, const double *__y, const double *__x, const double *__p, int give_log, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, int __ncovars, const double *__covars, double t)\n{\n",
  step.fn="\nvoid __pomp_stepfn (double *__x, const double *__p, const int *__stateindex, const int *__parindex, const int *__covindex, int __covdim, const double *__covars, double t, double dt)\n{\n",
  rate.fn="\ndouble __pomp_ratefn (int j, double t, double *__x, const double *__p, const int *__stateindex, const int *__parindex, const int *__covindex, int __covdim, const double *__covars)\n{\n  double rate = 0.0;  \n",
  dens.fn="\nvoid __pomp_densfn (double *__loglik, const double *__x1, const double *__x2, double t1, double t2, const double *__p, const int *__stateindex, const int *__parindex, const int *__covindex, int __ncovars, const double *__covars)\n{\n",
  skeleton="\nvoid __pomp_skelfn (double *__f, const double *__x, const double *__p, const int *__stateindex, const int *__parindex, const int *__covindex, int __ncovars, const double *__covars, double t)\n{\n",
  initializer="\nvoid __pomp_initializer (double *__x, const double *__p, double t, const int *__stateindex, const int *__parindex, const int *__covindex, const double *__covars)\n{\n",
  fromEstimationScale="\nvoid __pomp_par_trans (double *__pt, const double *__p, const int *__parindex)\n{\n",
  toEstimationScale="\nvoid __pomp_par_untrans (double *__pt, const double *__p, const int *__parindex)\n{\n",
  rprior="\nvoid __pomp_rprior (double *__p, const int *__parindex)\n{\n",
  dprior="\nvoid __pomp_dprior (double *__lik, const double *__p, int give_log, const int *__parindex)\n{\n"
)

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

registration <- "R_RegisterCCallable(\"{%name%}\", \"__pomp_{%fun%}\", (DL_FUNC) __pomp_{%fun%});\n"

decl <- list(
  periodic_bspline_basis_eval="\tvoid (*periodic_bspline_basis_eval)(double,double,int,int,double*);\nperiodic_bspline_basis_eval = (void (*)(double,double,int,int,double*)) R_GetCCallable(\"pomp\",\"periodic_bspline_basis_eval\");\n",
  get_pomp_userdata_int="\tconst int * (*get_pomp_userdata_int)(const char *);\nget_pomp_userdata_int = (const int *(*)(const char*)) R_GetCCallable(\"pomp\",\"get_pomp_userdata_int\");\n",
  get_pomp_userdata_double="\tconst double * (*get_pomp_userdata_double)(const char *);\nget_pomp_userdata_double = (const double *(*)(const char*)) R_GetCCallable(\"pomp\",\"get_pomp_userdata_double\");\n",
  `get_pomp_userdata(\\b|[^_])`="\tconst SEXP (*get_pomp_userdata)(const char *);\nget_pomp_userdata = (const SEXP (*)(const char*)) R_GetCCallable(\"pomp\",\"get_pomp_userdata\");\n"
)

footer <- list(
  rmeasure="\n}\n\n",
  dmeasure="\n}\n\n",
  step.fn="\n}\n\n",
  rate.fn="  return rate;\n}\n\n",
  dens.fn="\n}\n\n",
  skeleton="\n}\n\n",
  initializer="\n}\n\n",
  fromEstimationScale="\n}\n\n",
  toEstimationScale="\n}\n\n",
  rprior="\n}\n\n",
  dprior="\n}\n\n",
  globals="\n\n",
  registration="}\n\n"
)

stackhandling <- "
static int __pomp_load_stack = 0;\n
void __pomp_load_stack_incr (void) {
  ++__pomp_load_stack;
}\n
void __pomp_load_stack_decr (int *val) {
  *val = --__pomp_load_stack;
}\n"

## utility.fns <- list()
