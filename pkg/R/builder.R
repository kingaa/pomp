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
                         statenames, paramnames, zeronames,
                         rmeasure, dmeasure, step.fn, step.fn.delta.t,
                         skeleton, skeleton.type, skelmap.delta.t = 1,
                         link = TRUE) {
  obsnames <- names(data)
  obsnames <- setdiff(obsnames,times)
  solib <- pompCBuilder(
                        name=name,
                        statenames=statenames,
                        paramnames=paramnames,
                        obsnames=obsnames,
                        rmeasure=rmeasure,
                        dmeasure=dmeasure,
                        step.fn=step.fn,
                        skeleton=skeleton
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
       PACKAGE=name,
       obsnames=obsnames,
       statenames=statenames,
       paramnames=paramnames,
       zeronames=zeronames
       )
}

pompLink <- function (name) {
  solib <- paste(name,.Platform$dynlib.ext,sep="")
  dyn.load(solib)
}

pompUnlink <- function (name) {
  solib <- paste(name,.Platform$dynlib.ext,sep="")
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
               file="/* pomp model file: {%name%} */\n\n#include <pomp.h>\n#include <R_ext/Rdynload.h>\n",
               rmeasure="\nvoid {%name%}_rmeasure (double *__y, double *__x, double *__p, int *__obsindex, int *__stateindex, int *__parindex, int *__covindex, int __ncovars, double *__covars, double t)\n{\n",
               dmeasure= "\nvoid {%name%}_dmeasure (double *__lik, double *__y, double *__x, double *__p, int __give_log, int *__obsindex, int *__stateindex, int *__parindex, int *__covindex, int __ncovars, double *__covars, double t)\n{\n",
               step.fn="\nvoid {%name%}_stepfn (double *__x, const double *__p, const int *__stateindex, const int *__parindex, const int *__covindex, int __covdim, const double *__covar, double t, double dt)\n{\n",
               skeleton="\nvoid {%name%}_skelfn (double *__f, double *__x, double *__p, int *__stateindex, int *__parindex, int *__covindex, int __ncovars, double *__covars, double t)\n{\n"
               )

decl <- list(
             periodic_bspline_basis_eval="void (*periodic_bspline_basis_eval)(double,double,int,int,double*);\nperiodic_bspline_basis_eval = (void (*)(double,double,int,int,double*)) R_GetCCallable(\"pomp\",\"periodic_bspline_basis_eval\");\n",
             reulermultinom="void (*reulermultinom)(int,double,double*,double,double*);\nreulermultinom = (void (*)(int,double,double*,double,double*)) R_GetCCallable(\"pomp\",\"reulermultinom\");\n"
              )

footer <- list(
               rmeasure="\n}\n\n",
               dmeasure="\nlik = (__give_log) ? log(lik) : lik;\n}\n\n",
               step.fn="\n}\n\n",
               skeleton="\n}\n\n"
               )

pompCBuilder <- function (name, statenames, paramnames, obsnames, rmeasure, dmeasure, step.fn, skeleton) {

  name <- cleanForC(name)
  statenames <- cleanForC(statenames)
  paramnames <- cleanForC(paramnames)
  obsnames <- cleanForC(obsnames)

  modelfile <- paste(name,".c",sep="")
  solib <- paste(name,.Platform$dynlib.ext,sep="")

  out <- file(description=modelfile,open="w")
  
  cat(file=out,render(header$file,name=name))
  ## variable/parameter/observations definitions
  for (v in seq_along(paramnames)) {
    cat(file=out,render(define$var,variable=paramnames[v],ptr='__p',ilist='__parindex',index=v-1))
  }
  for (v in seq_along(statenames)) {
    cat(file=out,render(define$var,variable=statenames[v],ptr='__x',ilist='__stateindex',index=v-1))
  }
  for (v in seq_along(obsnames)) {
    cat(file=out,render(define$var,variable=obsnames[v],ptr='__y',ilist='__obsindex',index=v-1))
  }
  for (v in seq_along(statenames)) {
    cat(file=out,render(define$var,variable=paste("D",statenames[v],sep=""),ptr='__f',ilist='__stateindex',index=v-1))
  }
  cat(file=out,render(define$var.alt,variable="lik",ptr='__lik',index=0))
  
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
  for (v in seq_along(statenames)) {
    cat(file=out,render(undefine$var,variable=statenames[v]))
  }
  for (v in seq_along(obsnames)) {
    cat(file=out,render(undefine$var,variable=obsnames[v]))
  }
  for (v in seq_along(statenames)) {
    cat(file=out,render(undefine$var,variable=paste("D",statenames[v],sep="")))
  }
  close(out)

  cflags <- paste("PKG_CFLAGS=\"",
                  Sys.getenv("PKG_CFLAGS"),
                  " -I",system.file("include",package="pomp"),"\"",
                  sep="")

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
