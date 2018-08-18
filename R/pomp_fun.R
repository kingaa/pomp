##' The "pomp_fun" class
##'
##' Definition and methods of the \sQuote{pomp_fun} class
##'
##' The \sQuote{pomp_fun} class helps to settle common issues associated with
##' user-defined functions which can be defined either via R code or by a
##' native, compiled routine.  It is not exported to userland.
##'
##' @name pomp_fun
##' @rdname pomp_fun
##' @include pomp-package.R csnippet.R
##' @docType methods
##' @keywords programming internal
##'
##' @param f A function or the name of a native routine.
##' @param PACKAGE optional; the name of the dynamically-loadable library in
##' which the native function \code{f} can be found.
##' @param proto optional string; a prototype against which \code{f} will be
##' checked.
##' @param object,x the \sQuote{pomp_fun} object.
##' @author Aaron A. King
##' @seealso \code{\link{pomp}}
NULL

## also defined in 'pomp_defines.h'
pompfunmode <- list(undef=-1L,Rfun=0L,native=1L,regNative=2L)

setClass(
  "pomp_fun",
  slots=c(
    R.fun = "function",
    native.fun = "character",
    PACKAGE = "character",
    mode = "integer",
    address = "externalptr",
    obsnames = "character",
    statenames = "character",
    paramnames = "character",
    covarnames = "character",
    purpose = "character"
  ),
  prototype=prototype(
    R.fun=function (...) {
      stop("in ",sQuote("pomp_fun"),
        ": unreachable error: please report this bug!",call.=FALSE)
    },
    native.fun=character(0),
    PACKAGE=character(0),
    mode=pompfunmode$undef,
    obsnames = character(0),
    statenames = character(0),
    paramnames = character(0),
    covarnames = character(0),
    purpose = "a needed function"
  )
)

setGeneric(
  "pomp_fun",
  function (f, ...)
    standardGeneric("pomp_fun")
)

setMethod(
  "pomp_fun",
  signature=signature(f="missing"),
  definition=function (f, slotname = NULL, ...) {
    new("pomp_fun",purpose=as.character(slotname))
  }
)

setMethod(
  "pomp_fun",
  signature=signature(f="NULL"),
  definition=function (f, slotname = NULL, ...) {
    new("pomp_fun",purpose=as.character(slotname))
  }
)

setMethod(
  "pomp_fun",
  signature=signature(f="ANY"),
  definition=function (f, slotname = NULL, ...) {
    stop("bad option for ",sQuote(slotname)," argument",call.=FALSE)
  }
)

setMethod(
  "pomp_fun",
  signature=signature(f="function"),
  definition=function (f, proto = NULL, slotname = NULL, ...) {
    if (!is.null(proto)) {
      prototype <- as.character(proto)
      fname <- prototype[1]
      args <- prototype[-1]
      if (is.function(f)&&(!all(args%in%names(formals(f)))))
        stop("in ",sQuote(slotname),": ",
          sQuote(fname)," must be a function of prototype ",
          sQuote(deparse(proto)),call.=FALSE)
    }
    new("pomp_fun",R.fun=f,mode=pompfunmode$Rfun,purpose=as.character(slotname))
  }
)

setMethod(
  "pomp_fun",
  signature=signature(f="character"),
  definition=function (f, PACKAGE = NULL,
    obsnames = character(0), statenames = character(0),
    paramnames = character(0), covarnames = character(0),
    slotname = NULL, ...) {
    new(
      "pomp_fun",
      native.fun=f,
      PACKAGE=as.character(PACKAGE),
      mode=pompfunmode$native,
      obsnames=obsnames,
      statenames=statenames,
      paramnames=paramnames,
      covarnames=covarnames,
      purpose=as.character(slotname)
    )
  }
)

setMethod(
  "pomp_fun",
  signature=signature(f="Csnippet"),
  definition=function (f, slotname = NULL, libname = NULL,
    obsnames = character(0), statenames = character(0),
    paramnames = character(0), covarnames = character(0),
    Cname, ...) {
    slotname <- as.character(slotname)
    libname <- as.character(libname)
    new(
      "pomp_fun",
      native.fun=render(Cname,name=libname),
      PACKAGE=libname,
      mode=pompfunmode$regNative,
      obsnames=obsnames,
      statenames=statenames,
      paramnames=paramnames,
      covarnames=covarnames,
      purpose=slotname
    )
  }
)

setMethod(
  "pomp_fun",
  signature=signature(f="pomp_fun"),
  definition=function (f, ...) f
)
