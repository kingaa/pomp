##' The "pomp_fun" class
##'
##' Definition and methods of the \sQuote{pomp_fun} class.
##'
##' The \sQuote{pomp_fun} class implements a common interface for user-defined procedures that can be defined in terms of R code or by compiled native routines.
##'
##' @name pomp_fun
##' @rdname pomp_fun
##' @include pomp-package.R csnippet.R pstop.R undefined.R
##' @docType methods
##' @keywords programming internal
##' @concept extending the \pkg{pomp} package
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
pompfunmode <- list(undef=0L,Rfun=1L,native=2L,regNative=3L)

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
      pStop("pomp_fun","unreachable error: please report this bug!")
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

##' @rdname pomp_fun
setMethod(
  "pomp_fun",
  signature=signature(f="missing"),
  definition=function (slotname = NULL,
    obsnames = character(0), statenames = character(0),
    paramnames = character(0), covarnames = character(0), ...) {
    new("pomp_fun",
      obsnames=obsnames,
      statenames=statenames,
      paramnames=paramnames,
      covarnames=covarnames,
      purpose=as.character(slotname))
  }
)

setMethod(
  "pomp_fun",
  signature=signature(f="NULL"),
  definition=function (f, slotname = NULL,
    obsnames = character(0), statenames = character(0),
    paramnames = character(0), covarnames = character(0), ...) {
    new("pomp_fun",
      obsnames=obsnames,
      statenames=statenames,
      paramnames=paramnames,
      covarnames=covarnames,
      purpose=as.character(slotname))
  }
)

setMethod(
  "pomp_fun",
  signature=signature(f="ANY"),
  definition=function (f, slotname = NULL, ...) {
    pStop_("bad option for ",sQuote(slotname)," argument.")
  }
)

##' @rdname pomp_fun
setMethod(
  "pomp_fun",
  signature=signature(f="function"),
  definition=function (f, proto = NULL, slotname = NULL, ...) {
    if (!is.null(proto)) {
      prototype <- as.character(proto)
      fname <- prototype[1]
      args <- prototype[-1]
      if (is.function(f)&&(!all(args%in%names(formals(f)))))
        pStop(slotname,
          sQuote(fname)," must be a function of prototype ",
          sQuote(deparse(proto)))
    }
    new("pomp_fun",R.fun=f,mode=pompfunmode$Rfun,purpose=as.character(slotname))
  }
)

##' @rdname pomp_fun
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

##' @rdname pomp_fun
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

##' @rdname pomp_fun
setMethod(
  "pomp_fun",
  signature=signature(f="pomp_fun"),
  definition=function (f, ...) f
)

setMethod(
  "show",
  signature=signature("pomp_fun"),
  definition=function (object) {

    mode <- object@mode

    if (mode==pompfunmode$Rfun) { # R function

      cat("R function\n  - definition: ")
      f <- object@R.fun
      environment(f) <- globalenv()
      print(f)

    } else if (mode==pompfunmode$native) { # user supplied native code

      cat("native function\n  - name: ",sQuote(object@native.fun),"\n",sep="")
      if (length(object@PACKAGE)>0)
        cat("  - dynamically loaded from: ",
          sQuote(object@PACKAGE),sep="")
      cat("\n")

    } else if (mode==pompfunmode$regNative) { # built from C snippets

      cat("native function\n  - name: ",sQuote(object@native.fun),"\n",sep="")
      if (length(object@PACKAGE)>0)
        cat("  - defined by a C snippet in library ",
          sQuote(object@PACKAGE),sep="")
      cat("\n")

    } else {

      cat("<default>")

    }
    cat("\n")
  }
)

##' @rdname undefined
setMethod(
  "undefined",
  signature=signature(object="pomp_fun"),
  definition=function (object, ...) {
    object@mode == pompfunmode$undef
  }
)
