setClass(
         "pompCode",
         representation=representation(
           type="character",
           slot="character",
           text="character",
           fun="function"
           ),
         prototype=prototype(
           type="ccode",
           slot=character(0),
           text=character(0),
           fun=function(...)stop("function not specified")
           )
         )


CCode <- function (text, slot) {
  new("pompCode",type="ccode",slot=as.character(slot))
}

pompBuilder <- function (data, times, t0, name,
                         statenames, paramnames, tcovar, covar,
                         rmeasure, dmeasure, step.fn, step.fn.delta.t,
                         skeleton, skeleton.type, skelmap.delta.t = 1,
                         parameter.transform, parameter.inv.transform,
                         ..., link = TRUE) {
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
  solib <- pompCBuilder(
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
                        parameter.inv.transform=parameter.inv.transform
                        )
  if (link) pompLink(name)
  pomp(
       data=data,times=times,t0=t0,
       rprocess=euler.sim(
         step.fun=render("{%name%}_stepfn",name=name),
         delta.t=step.fn.delta.t,
         PACKAGE=name
         ),
       rmeasure=render("{%name%}_rmeasure",name=name),
       dmeasure=render("{%name%}_dmeasure",name=name),
       skeleton=render("{%name%}_skelfn",name=name),
       skeleton.type=skeleton.type,
       skelmap.delta.t=skelmap.delta.t,
       parameter.transform=render("{%name%}_par_trans",name=name),
       parameter.inv.transform=render("{%name%}_par_untrans",name=name),
       PACKAGE=name,
       obsnames=obsnames,
       statenames=statenames,
       paramnames=paramnames,
       covarnames=covarnames,
       tcovar=tcovar,
       covar=covar,
       ...
       )
}

pompLink <- function (name) {
  solib <- paste0(name,.Platform$dynlib.ext)
  dyn.load(solib)
}

pompUnlink <- function (name) {
  solib <- paste0(name,.Platform$dynlib.ext)
  dyn.unload(solib)
}

define <- list(
               var="#define {%variable%}\t({%ptr%}[{%ilist%}[{%index%}]])\n",
               var.alt="#define {%variable%}\t({%ptr%}[{%index%}])\n"
               )

undefine <- list(
                 var="#undef {%variable%}\n"
                 )

header <- list(
               file="/* pomp model file: {%name%} */\n\n#include <pomp.h>\n#include <R_ext/Rdynload.h>\n\n",
               rmeasure="\nvoid {%name%}_rmeasure (double *__y, double *__x, double *__p, int *__obsindex, int *__stateindex, int *__parindex, int *__covindex, int __ncovars, double *__covars, double t)\n{\n",
               dmeasure= "\nvoid {%name%}_dmeasure (double *__lik, double *__y, double *__x, double *__p, int give_log, int *__obsindex, int *__stateindex, int *__parindex, int *__covindex, int __ncovars, double *__covars, double t)\n{\n",
               step.fn="\nvoid {%name%}_stepfn (double *__x, const double *__p, const int *__stateindex, const int *__parindex, const int *__covindex, int __covdim, const double *__covars, double t, double dt)\n{\n",
               skeleton="\nvoid {%name%}_skelfn (double *__f, double *__x, double *__p, int *__stateindex, int *__parindex, int *__covindex, int __ncovars, double *__covars, double t)\n{\n",
               parameter.transform="\nvoid {%name%}_par_trans (double *__pt, double *__p, int *__parindex)\n{\n",
               parameter.inv.transform="\nvoid {%name%}_par_untrans (double *__pt, double *__p, int *__parindex)\n{\n"
               )

decl <- list(
             periodic_bspline_basis_eval="\tvoid (*periodic_bspline_basis_eval)(double,double,int,int,double*);\nperiodic_bspline_basis_eval = (void (*)(double,double,int,int,double*)) R_GetCCallable(\"pomp\",\"periodic_bspline_basis_eval\");\n",
             reulermultinom="\tvoid (*reulermultinom)(int,double,double*,double,double*);\nreulermultinom = (void (*)(int,double,double*,double,double*)) R_GetCCallable(\"pomp\",\"reulermultinom\");\n",
             get_pomp_userdata="\tconst SEXP (*get_pomp_userdata)(const char *);\npomp_get_userdata = (const SEXP (*)(const char*)) R_GetCCallable(\"pomp\",\"get_pomp_userdata\");\n",
             get_pomp_userdata_int="\tconst int * (*get_pomp_userdata_int)(const char *);\npomp_get_userdata_int = (const int *(*)(const char*)) R_GetCCallable(\"pomp\",\"get_pomp_userdata_int\");\n",
             get_pomp_userdata_double="\tconst double * (*get_pomp_userdata_double)(const char *);\npomp_get_userdata_double = (const double *(*)(const char*)) R_GetCCallable(\"pomp\",\"get_pomp_userdata_double\");\n"
             )

footer <- list(
               rmeasure="\n}\n\n",
               dmeasure="\n}\n\n",
               step.fn="\n}\n\n",
               skeleton="\n}\n\n",
               parameter.transform="\n}\n\n",
               parameter.inv.transform="\n}\n\n"
               )

utility.fns <- list(
                    )


pompCBuilder <- function (name, statenames, paramnames, covarnames, obsnames, rmeasure, dmeasure,
                          step.fn, skeleton, parameter.transform, parameter.inv.transform)
{
  if (missing(name)) stop(sQuote("name")," must be supplied");
  if (missing(statenames)) stop(sQuote("statenames")," must be supplied");
  if (missing(paramnames)) stop(sQuote("paramnames")," must be supplied");
  if (missing(obsnames)) stop(sQuote("obsnames")," must be supplied");
  if (missing(covarnames)) stop(sQuote("covarnames")," must be supplied");

  mpt <- missing(parameter.transform)
  mpit <- missing(parameter.inv.transform)
  if (xor(mpt,mpit))
    stop("if you supply one transformation function, you must supply its inverse")
  has.trans <- !mpt

  name <- cleanForC(name)
  statenames <- cleanForC(statenames)
  paramnames <- cleanForC(paramnames)
  covarnames <- cleanForC(covarnames)
  obsnames <- cleanForC(obsnames)

  modelfile <- paste0(name,".c")
  solib <- paste0(name,.Platform$dynlib.ext)

  out <- file(description=modelfile,open="w")
  
  cat(file=out,render(header$file,name=name))

  for (f in utility.fns) {
    cat(file=out,f)
  }

  ## variable/parameter/observations definitions
  for (v in seq_along(paramnames)) {
    cat(file=out,render(define$var,variable=paramnames[v],ptr='__p',ilist='__parindex',index=v-1))
  }
  for (v in seq_along(covarnames)) {
    cat(file=out,render(define$var,variable=covarnames[v],ptr='__covars',ilist='__covindex',index=v-1))
  }
  for (v in seq_along(statenames)) {
    cat(file=out,render(define$var,variable=statenames[v],ptr='__x',ilist='__stateindex',index=v-1))
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

  if (has.trans) {
    ## parameter transformation function
    cat(file=out,render(header$parameter.transform,name=name))
    for (fn in names(decl)) {
      if (grepl(fn,parameter.transform))
        cat(file=out,decl[[fn]])
    }
    cat(file=out,parameter.transform,footer$parameter.transform)
    ## inverse parameter transformation function
    cat(file=out,render(header$parameter.inv.transform,name=name))
    for (fn in names(decl)) {
      if (grepl(fn,parameter.inv.transform))
        cat(file=out,decl[[fn]])
    }
    cat(file=out,parameter.inv.transform,footer$parameter.inv.transform)
  }

  ## rmeasure function
  cat(file=out,render(header$rmeasure,name=name),rmeasure,footer$rmeasure)

  ## dmeasure function
  cat(file=out,render(header$dmeasure,name=name),dmeasure,footer$dmeasure)

  ## Euler step function
  cat(file=out,render(header$step.fn,name=name))
  for (fn in names(decl)) {
    if (grepl(fn,step.fn))
      cat(file=out,decl[[fn]])
  }
  cat(file=out,step.fn,footer$step.fn)

  ## skeleton function
  cat(file=out,render(header$skeleton,name=name))
  for (fn in names(decl)) {
    if (grepl(fn,skeleton))
      cat(file=out,decl[[fn]])
  }
  cat(file=out,skeleton,footer$skeleton)

  ## undefine variables
  for (v in seq_along(paramnames)) {
    cat(file=out,render(undefine$var,variable=paramnames[v]))
  }
  for (v in seq_along(covarnames)) {
    cat(file=out,render(undefine$var,variable=covarnames[v]))
  }
  for (v in seq_along(statenames)) {
    cat(file=out,render(undefine$var,variable=statenames[v]))
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
  close(out)

  cflags <- paste0("PKG_CFLAGS=\"",
                  Sys.getenv("PKG_CFLAGS"),
                  " -I",system.file("include",package="pomp"),"\"")

  rv <- system2(
                command=R.home("bin/R"),
                args=c("CMD","SHLIB","-o",solib,modelfile),
                env=cflags
                )
  if (rv!=0)
    stop("cannot compile shared-object library ",sQuote(solib))
  else
    cat("model codes written to",sQuote(modelfile),"\nlink to shared-object library",sQuote(solib),"\n")

  invisible(solib)
}

cleanForC <- function (text) {
  text <- as.character(text)
  text <- gsub("\\.","_",text)
  text
}

render <- function (template, ...) {
  vars=list(...)
  n <- sapply(vars,length)
  if (!all((n==max(n))|(n==1))) stop("incommensurate lengths of replacements")
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
