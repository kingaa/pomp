pompCBuilder <- function (name = NULL, dir = NULL,
                          statenames, paramnames, covarnames, obsnames,
                          rmeasure, dmeasure, step.fn, skeleton,
                          parameter.transform, parameter.inv.transform,
                          rprior, dprior, globals,
                          verbose = getOption("verbose",FALSE))
{

  if (is.null(name)) name <- randomName()
  if (is.null(dir)) dir <- tempdir()

  if (missing(globals)) globals <- character(0)

  name <- cleanForC(name)
  statenames <- cleanForC(statenames)
  paramnames <- cleanForC(paramnames)
  covarnames <- cleanForC(covarnames)
  obsnames <- cleanForC(obsnames)

  stem <- file.path(dir,name)
  if (.Platform$OS.type=="windows") {
    stem <- gsub("\\","/",stem,fixed=TRUE)
  }
  modelfile <- paste0(stem,".c") 
  solib <- paste0(stem,.Platform$dynlib.ext)

  if (.Platform$OS.type=="unix") {
    pompheader <- "pomp.h"
  } else {
    pompheader <- system.file("include/pomp.h",package="pomp")
  }
  
  out <- file(description=modelfile,open="w")
  
  cat(file=out,render(header$file,name=name,pompheader=pompheader))

  for (f in utility.fns) {
    cat(file=out,f)
  }

  cat(file=out,globals,footer$globals)

  ## variable/parameter/observations definitions
  for (v in seq_along(paramnames)) {
    cat(file=out,render(define$var,variable=paramnames[v],ptr='__p',ilist='__parindex',index=v-1))
  }
  for (v in seq_along(statenames)) {
    cat(file=out,render(define$var,variable=statenames[v],ptr='__x',ilist='__stateindex',index=v-1))
  }
  for (v in seq_along(covarnames)) {
    cat(file=out,render(define$var,variable=covarnames[v],ptr='__covars',ilist='__covindex',index=v-1))
  }
  for (v in seq_along(obsnames)) {
    cat(file=out,render(define$var,variable=obsnames[v],ptr='__y',ilist='__obsindex',index=v-1))
  }
  for (v in seq_along(statenames)) {
    cat(file=out,render(define$var,variable=paste0("D",statenames[v]),ptr='__f',ilist='__stateindex',index=v-1))
  }
  for (v in seq_along(paramnames)) {
    cat(file=out,render(define$var,variable=paste0("T",paramnames[v]),ptr='__pt',ilist='__parindex',index=v-1))
  }
  cat(file=out,render(define$var.alt,variable="lik",ptr='__lik',index=0))

  ## list of functions to register
  registry <- c("load_stack_incr","load_stack_decr")

  ## parameter transformation function
  if (!missing(parameter.transform)) {
    registry <- c(registry,"par_trans")
    cat(file=out,render(header$parameter.transform,name=name))
    cat(file=out,callable.decl(parameter.transform))
    cat(file=out,parameter.transform,footer$parameter.transform)
  }

  ## inverse parameter transformation function
  if (!missing(parameter.inv.transform)) {
    registry <- c(registry,"par_untrans")
    cat(file=out,render(header$parameter.inv.transform,name=name))
    cat(file=out,callable.decl(parameter.inv.transform))
    cat(file=out,parameter.inv.transform,footer$parameter.inv.transform)
  }

  ## rmeasure function
  if (!missing(rmeasure)) {
    registry <- c(registry,"rmeasure")
    cat(file=out,render(header$rmeasure,name=name),rmeasure,footer$rmeasure)
  }

  ## dmeasure function
  if (!missing(dmeasure)) {
    registry <- c(registry,"dmeasure")
    cat(file=out,render(header$dmeasure,name=name),dmeasure,footer$dmeasure)
  }

  ## Euler step function
  if (!missing(step.fn)) {
    registry <- c(registry,"stepfn")
    cat(file=out,render(header$step.fn,name=name))
    cat(file=out,callable.decl(step.fn))
    cat(file=out,step.fn,footer$step.fn)
  }

  ## skeleton function
  if (!missing(skeleton)) {
    registry <- c(registry,"skelfn")
    cat(file=out,render(header$skeleton,name=name))
    cat(file=out,callable.decl(skeleton))
    cat(file=out,skeleton,footer$skeleton)
  }

  ## rprior function
  if (!missing(rprior)) {
    registry <- c(registry,"rprior")
    cat(file=out,render(header$rprior,name=name),rprior,footer$rprior)
  }

  ## dprior function
  if (!missing(dprior)) {
    registry <- c(registry,"dprior")
    cat(file=out,render(header$dprior,name=name),dprior,footer$dprior)
  }

  ## undefine variables
  for (v in seq_along(paramnames)) {
    cat(file=out,render(undefine$var,variable=paramnames[v]))
  }
  for (v in seq_along(statenames)) {
    cat(file=out,render(undefine$var,variable=statenames[v]))
  }
  for (v in seq_along(covarnames)) {
    cat(file=out,render(undefine$var,variable=covarnames[v]))
  }
  for (v in seq_along(obsnames)) {
    cat(file=out,render(undefine$var,variable=obsnames[v]))
  }
  for (v in seq_along(statenames)) {
    cat(file=out,render(undefine$var,variable=paste0("D",statenames[v])))
  }
  for (v in seq_along(paramnames)) {
    cat(file=out,render(undefine$var,variable=paste0("T",paramnames[v])))
  }

  cat(file=out,stackhandling)

  ## registration
  cat(file=out,render(header$registration,name=name))
  for (v in registry)
    cat(file=out,render(registration,name=name,fun=v))
  cat(file=out,footer$registration)

  close(out)

  cflags <- paste0("PKG_CFLAGS=\"",
                  Sys.getenv("PKG_CFLAGS"),
                  " -I",system.file("include",package="pomp"),"\"")

  rv <- system2(
                command=R.home("bin/R"),
                args=c("CMD","SHLIB","-o",solib,modelfile),
                env=cflags,
                stdout=if (verbose) "" else NULL
                )
  if (rv!=0)
    stop("cannot compile shared-object library ",sQuote(solib))
  else if (verbose)
    cat("model codes written to",sQuote(modelfile),
        "\nlink to shared-object library",sQuote(solib),"\n")

  invisible(c(name,solib))
}

callable.decl <- function (code) {
  fns <- vapply(names(decl),grepl,logical(1),code,perl=TRUE)
  do.call(paste0,decl[fns])
}

missing.fun <- function (name) {
  paste0("  error(\"'",name,"' not defined\");")
}

randomName <- function (stem = "pomp", size = 2) {
  paste0(stem,
         paste(
               format(
                      as.hexmode(ceiling(runif(n=size,max=2^24))),
                      upper.case=TRUE
                      ),
               collapse=""
               )
         )
}

cleanForC <- function (text) {
  text <- as.character(text)
  text <- gsub("\\.","_",text)
  text
}

render <- function (template, ...) {
  vars=list(...)
  n <- sapply(vars,length)
  if (!all((n==max(n))|(n==1)))
    stop("incommensurate lengths of replacements")
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

pompBuilder <- function (data, times, t0, name,
                         statenames, paramnames, tcovar, covar,
                         rmeasure, dmeasure, step.fn, step.fn.delta.t,
                         skeleton, skeleton.type = c("map","vectorfield"),
                         skelmap.delta.t = 1,
                         parameter.transform, parameter.inv.transform,
                         rprior, dprior,
                         globals, ..., save = FALSE) {
  
  if (!is.data.frame(data)) stop(sQuote("data")," must be a data-frame")
  obsnames <- names(data)
  obsnames <- setdiff(obsnames,times)
  if (!missing(covar)) {
    if (!is.data.frame(covar)) stop(sQuote("covar")," must be a data-frame")
    covarnames <- colnames(covar)
    covarnames <- setdiff(covarnames,tcovar)
  } else {
    covar <- matrix(data=0,nrow=0,ncol=0)
    tcovar <- numeric(0)
    covarnames <- character(0)
  }
  skeleton.type <- match.arg(skeleton.type)

  if (missing(statenames)) stop(sQuote("statenames")," must be supplied");
  if (missing(paramnames)) stop(sQuote("paramnames")," must be supplied");

  mpt <- missing(parameter.transform)
  mpit <- missing(parameter.inv.transform)
  if (xor(mpt,mpit))
    stop("if you supply one transformation function, you must supply its inverse")

  pompCBuilder(
               name=name,
               statenames=statenames,
               paramnames=paramnames,
               covarnames=covarnames,
               obsnames=obsnames,
               rmeasure=rmeasure,
               dmeasure=dmeasure,
               step.fn=step.fn,
               skeleton=skeleton,
               parameter.transform=parameter.transform,
               parameter.inv.transform=parameter.inv.transform,
               rprior=rprior,
               dprior=dprior,
               globals=globals,
               dir=if (save) getwd() else NULL
               ) -> solib

  name <- solib[1]

  pomp(
       data=data,
       times=times,
       t0=t0,
       rprocess=euler.sim(
         step.fun=render(fnames$step.fn,name=name),
         delta.t=step.fn.delta.t,
         PACKAGE=name
         ),
       rmeasure=render(fnames$rmeasure,name=name),
       dmeasure=render(fnames$dmeasure,name=name),
       skeleton=render(fnames$skeleton,name=name),
       skeleton.type=skeleton.type,
       skelmap.delta.t=skelmap.delta.t,
       parameter.transform=render(fnames$parameter.transform,name=name),
       parameter.inv.transform=render(fnames$parameter.inv.transform,name=name),
       rprior=render(fnames$rprior,name=name),
       dprior=render(fnames$dprior,name=name),
       PACKAGE=name,
       statenames=statenames,
       paramnames=paramnames,
       tcovar=tcovar,
       covar=covar,
       ...,
       .solibfile=list(solib)
       )
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
               file="/* pomp model file: {%name%} */\n\n#include <{%pompheader%}>\n#include <R_ext/Rdynload.h>\n\n",
               registration="\nvoid R_init_{%name%} (DllInfo *info)\n{\n",
               rmeasure="\nvoid __pomp_rmeasure (double *__y, double *__x, double *__p, int *__obsindex, int *__stateindex, int *__parindex, int *__covindex, int __ncovars, double *__covars, double t)\n{\n",
               dmeasure= "\nvoid __pomp_dmeasure (double *__lik, double *__y, double *__x, double *__p, int give_log, int *__obsindex, int *__stateindex, int *__parindex, int *__covindex, int __ncovars, double *__covars, double t)\n{\n",
               step.fn="\nvoid __pomp_stepfn (double *__x, const double *__p, const int *__stateindex, const int *__parindex, const int *__covindex, int __covdim, const double *__covars, double t, double dt)\n{\n",
               skeleton="\nvoid __pomp_skelfn (double *__f, double *__x, double *__p, int *__stateindex, int *__parindex, int *__covindex, int __ncovars, double *__covars, double t)\n{\n",
               parameter.transform="\nvoid __pomp_par_trans (double *__pt, double *__p, int *__parindex)\n{\n",
               parameter.inv.transform="\nvoid __pomp_par_untrans (double *__pt, double *__p, int *__parindex)\n{\n",
               rprior="\nvoid __pomp_rprior (double *__p, int *__parindex)\n{\n",
               dprior="\nvoid __pomp_dprior (double *__lik, double *__p, int give_log, int *__parindex)\n{\n"
               )

fnames <- list(
               rmeasure="__pomp_rmeasure",
               dmeasure= "__pomp_dmeasure",
               step.fn="__pomp_stepfn",
               skeleton="__pomp_skelfn",
               parameter.transform="__pomp_par_trans",
               parameter.inv.transform="__pomp_par_untrans",
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
               skeleton="\n}\n\n",
               parameter.transform="\n}\n\n",
               parameter.inv.transform="\n}\n\n",
               rprior="\n}\n\n",
               dprior="\n}\n\n",
               globals="\n",
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

utility.fns <- list()
