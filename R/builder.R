##' Hitching C snippets and R functions to pomp_fun objects
##'
##' The algorithms in \pkg{pomp} are formulated using \R functions that access the \link[=basic_components]{basic model components}
##' (\code{rprocess}, \code{dprocess}, \code{rmeasure}, \code{dmeasure}, etc.).
##' For short, we refer to these elementary functions as \dQuote{\link{workhorses}}.
##' In implementing a model, the user specifies basic model components
##' using functions, procedures in dynamically-linked libraries, or C snippets.
##' Each component is then packaged into a \sQuote{pomp_fun} objects, which gives a uniform interface.
##' The construction of \sQuote{pomp_fun} objects is handled by the \code{hitch} function,
##' which conceptually \dQuote{hitches} the workhorses to the user-defined procedures.
##'
##' @name hitch
##' @docType methods
##' @include pomp_class.R csnippet.R safecall.R templates.R
##'
##' @importFrom digest digest
##' @importFrom stats runif
##'
##' @param \dots named arguments representing the user procedures to be hitched.
##' These can be functions, character strings naming routines in external, dynamically-linked libraries, C snippets, or \code{NULL}.
##' The first three are converted by \code{hitch} to \sQuote{pomp_fun} objects which perform the indicated computations.
##' \code{NULL} arguments are translated to default \sQuote{pomp_fun} objects.
##' If any of these procedures are already \sQuote{pomp_fun} objects, they are returned unchanged.
##'
##' @param templates named list of templates.
##' Each workhorse must have a corresponding template.
##' See \code{pomp:::workhorse_templates} for a list.
##'
##' @param obsnames,statenames,paramnames,covarnames character vectors
##' specifying the names of observable variables, latent state variables, parameters, and covariates, respectively.
##' These are only needed if one or more of the horses are furnished as C snippets.
##'
##' @param PACKAGE optional character;
##' the name (without extension) of the external, dynamically loaded library in which any native routines are to be found.
##' This is only useful if one or more of the model components has been specified using a precompiled dynamically loaded library;
##' it is not used for any component specified using C snippets.
##' \code{PACKAGE} can name at most one library.
##'
##' @param globals optional character;
##' arbitrary C code that will be hard-coded into the shared-object library created when C snippets are provided.
##' If no C snippets are used, \code{globals} has no effect.
##'
##' @param cdir optional character variable.
##' \code{cdir} specifies the name of the directory within which C snippet code will be compiled.
##' By default, this is in a temporary directory specific to the \R session.
##' One can also set this directory using the \code{pomp_cdir} global option.
##'
##' @param cfile optional character variable.
##' \code{cfile} gives the name of the file (in directory \code{cdir}) into which C snippet codes will be written.
##' By default, a random filename is used.
##' If the chosen filename would result in over-writing an existing file, an error is generated.
##'
##' @param shlib.args optional character variables.
##' Command-line arguments to the \code{R CMD SHLIB} call that compiles the C snippets.
##' One can, for example, specify libraries against which the C snippets are to be linked.
##' In doing so, take care to make sure the appropriate header files are available to the C snippets, e.g., using the \code{globals} argument.
##' See \code{\link{Csnippet}} for more information.
##'
##' @param compile logical;
##' if \code{FALSE}, compilation of the C snippets will be postponed until they are needed.
##'
##' @param verbose logical.
##' Setting \code{verbose=TRUE} will cause additional information to be displayed.
##'
##' @return
##' \code{hitch} returns a named list of length two.  The element named
##' \dQuote{funs} is itself a named list of \sQuote{pomp_fun} objects, each of
##' which corresponds to one of the horses passed in.  The element named
##' \dQuote{lib} contains information on the shared-object library created
##' using the C snippets (if any were passed to \code{hitch}).  If no C
##' snippets were passed to \code{hitch}, \code{lib} is \code{NULL}.
##' Otherwise, it is a length-3 named list with the following elements:
##' \describe{
##' \item{name}{The name of the library created.}
##' \item{dir}{ The
##' directory in which the library was created.  If this is \code{NULL}, the
##' library was created in the session's temporary directory.  }
##' \item{src}{ A
##' character string with the full contents of the C snippet file.  } }
##'
##' @author Aaron A. King
##'
##' @seealso \code{\link{pomp}}, \code{\link{spy}}
##' @concept extending the pomp package
##' @concept low-level interface
##'
NULL

## 'hitch' takes (as '...') the workhorse specifications (as R functions or
## C snippets, processes these, and hitches them to 'pomp_fun' objects
## suitable for use in the appropriate slots in 'pomp' objects.

##' @rdname hitch
##' @export
hitch <- function (..., templates,
  obsnames, statenames, paramnames, covarnames,
  PACKAGE, globals, cfile, cdir = getOption("pomp_cdir", NULL),
  shlib.args, compile = TRUE,
  verbose = getOption("verbose", FALSE)) {

  if (missing(templates))
    pStop(sQuote("templates")," must be supplied.")

  if (missing(cfile)) cfile <- NULL
  if (missing(cdir)) cdir <- NULL
  if (missing(PACKAGE)) PACKAGE <- NULL
  if (missing(globals)) globals <- NULL
  if (missing(shlib.args)) shlib.args <- NULL
  PACKAGE <- as.character(PACKAGE)
  compile <- as.logical(compile)

  ## defaults for names of states, parameters, observations, and covariates
  if (missing(statenames)) statenames <- NULL
  if (missing(paramnames)) paramnames <- NULL
  if (missing(obsnames)) obsnames <- NULL
  if (missing(covarnames)) covarnames <- NULL

  statenames <- as.character(statenames)
  paramnames <- as.character(paramnames)
  obsnames <- as.character(obsnames)
  covarnames <- as.character(covarnames)

  if (anyDuplicated(c(statenames,paramnames,obsnames,covarnames))) {
    pStop("the variable names in ",sQuote("statenames"),", ",
      sQuote("paramnames"),", ",sQuote("covarnames"),
      ", and ",sQuote("obsnames")," must be unique and disjoint.")
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
            compile=compile,
            verbose=verbose
          ),
          snippets
        )
      ),
      error = function (e)
        pStop_("error in building shared-object library from C snippets: ",conditionMessage(e)) #nocov
    )
    libname <- lib$name
  } else {
    libname <- ""
    lib <- NULL
  }

  funs <- vector(mode="list",length=length(horses))
  names(funs) <- names(horses)

  for (s in names(funs)) {
    funs[[s]] <- pomp_fun(
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
  globals, shlib.args = NULL, compile,
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
  utils <- which(
    vapply(
      seq_along(pomp_templates$utilities),
      \(x) any(grepl(pomp_templates$utilities[[x]]$trigger,
        snippets)),
      logical(1L)
    )
  )

  ## rely on "-I" flags under *nix
  if (.Platform$OS.type=="unix") {
    pompheader <- "pomp.h"
  } else {
    pompheader <- system.file("include/pomp.h",package="pomp") #nocov
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
      compile=compile,
      verbose=verbose
    ),
    error = function (e)
      pStop(who="Cbuilder","compilation error: ",conditionMessage(e)) #nocov
  )

  invisible(list(name=name,dir=dir,src=csrc))
}

pompCompile <- function (fname, direc, src, shlib.args = NULL,
  compile = TRUE, verbose) {

  stem <- file.path(direc,fname)
  if (.Platform$OS.type=="windows")
    stem <- gsub("\\","/",stem,fixed=TRUE) #nocov

  modelfile <- paste0(stem,".c")

  tryCatch(
    cat(src,file=modelfile),
    error = function (e) {
      pStop_("cannot write file ",sQuote(modelfile))   #nocov
    }
  )

  if (verbose) cat("model codes written to",sQuote(modelfile),"\n")

  cflags <- Sys.getenv("PKG_CPPFLAGS")
  cflags <- paste0("PKG_CPPFLAGS=\"",
    if (nchar(cflags)>0) paste0(cflags," ") else "",
    "-I",shQuote(system.file("include",package="pomp")),
    " -I",shQuote(getwd()),"\"")

  shlib.args <- as.character(shlib.args)

  solib <- paste0(stem,.Platform$dynlib.ext)

  if (compile) {

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
      error = function (e) pStop_("error compiling C snippets: ",conditionMessage(e)) #nocov
    )

    stat <- as.integer(attr(rv,"status"))

    if (length(stat) > 0 && stat != 0L) {
      pStop_("cannot compile shared-object library ",sQuote(solib),          #nocov
        ": status = ",stat,"\ncompiler messages:\n",paste(rv,collapse="\n")) #nocov
    } else if (verbose) {
      cat("compiler messages:",rv,sep="\n")
    }

  }

  invisible(solib)
}

srcDir <- function (dir, verbose) {
  if (is.null(dir)) {
    dir <- file.path(tempdir(),Sys.getpid())
  }
  if (!dir.exists(dir)) {
    if (verbose) cat("creating C snippet directory ",sQuote(dir),"\n") #nocov
    tryCatch(
      {
        dir.create(dir,recursive=TRUE,showWarnings=FALSE,mode="0700")
        stopifnot(dir.exists(dir))
      },
      error = function (e) pStop_("cannot create cache directory ",sQuote(dir))   #nocov
    )
  }
  dir
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
  n <- vapply(vars,length,integer(1L))
  if (!all((n==max(n))|(n==1)))
    pStop("incommensurate lengths of replacements.") #nocov
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
  define=r"{
#define {%variable%}  ({%cref%})
  }",
  undefine=r"{
#undef {%variable%}
  }",
  file=list(
    header=r"{
/* pomp C snippet file: {%name%} */
/* Time: {%timestamp%} */
/* Salt: {%salt%} */

#include <{%pompheader%}>
#include <R_ext/Rdynload.h>

    }"
    ),
  utilities=list(
    periodic_bspline_basis=list(
      trigger="periodic_bspline_basis_eval",
      header=r"{
static periodic_bspline_basis_eval_t *__pomp_periodic_bspline_basis_eval;
#define periodic_bspline_basis_eval(X,Y,M,N,Z)  (__pomp_periodic_bspline_basis_eval((X),(Y),(M),(N),(Z)))
      }",
      reg=r"{
__pomp_periodic_bspline_basis_eval = (periodic_bspline_basis_eval_t *) R_GetCCallable("pomp","periodic_bspline_basis_eval");
      }"
    ),
    get_userdata_int=list(
      trigger="get_userdata_int",
      header=r"{
static get_userdata_int_t *__pomp_get_userdata_int;
#define get_userdata_int(X)  (__pomp_get_userdata_int(X))
      }",
      reg=r"{
__pomp_get_userdata_int = (get_userdata_t *) R_GetCCallable("pomp","get_userdata_int");
      }"
    ),
    get_userdata_double=list(
      trigger="get_userdata_double",
      header=r"{
static get_userdata_double_t *__pomp_get_userdata_double;
#define get_userdata_double(X)  (__pomp_get_userdata_double(X))
      }",
      reg=r"{
__pomp_get_userdata_double = (get_userdata_double_t *) R_GetCCallable("pomp","get_userdata_double");
      }"
    ),
    get_userdata=list(
      trigger=r"{get_userdata(\b|[^_])}",
      header=r"{
static get_userdata_t *__pomp_get_userdata;
#define get_userdata(X)  (__pomp_get_userdata(X))
      }",
      reg=r"{
__pomp_get_userdata = (get_userdata_t *) R_GetCCallable("pomp","get_userdata");
      }"
    )
  ),
  stackhandling=r"{
static int __pomp_load_stack = 0;
void __pomp_load_stack_incr (void) {++__pomp_load_stack;}
void __pomp_load_stack_decr (int *val) {*val = --__pomp_load_stack;}
  }",
  registration=list(
    header=r"{
void R_init_{%name%} (DllInfo *info) {
    }",
    main=r"(
R_RegisterCCallable("{%name%}", "{%fun%}", (DL_FUNC) {%fun%});
    )",
    footer=r"{
}
    }"
  )
)
