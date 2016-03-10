## a class for functions that may be defined in R,
## using pre-written native routines, or C snippets

pompfunmode <- list(undef=-1L,Rfun=0L,native=1L,regNative=2L)

setClass(
         'pomp.fun',
         slots=c(
           R.fun = 'function',
           native.fun = 'character',
           PACKAGE = 'character',
           mode = 'integer',
           address = 'externalptr',
           obsnames = 'character',
           statenames = 'character',
           paramnames = 'character',
           covarnames = 'character'
           ),
         prototype=prototype(
           R.fun=function (...) {
             stop("unreachable error: please report this bug!")
           },
           native.fun=character(0),
           PACKAGE=character(0),
           mode=pompfunmode$undef,
           obsnames = character(0),
           statenames = character(0),
           paramnames = character(0),
           covarnames = character(0)
           )
         )

setMethod(
          "pomp.fun",
          signature=signature(f="missing"),
          definition=function (f, ...) {
            new("pomp.fun")
          }
          )

setMethod(
          "pomp.fun",
          signature=signature(f="NULL"),
          definition=function (f, ...) {
            new("pomp.fun")
          }
          )

setMethod(
          "pomp.fun",
          signature=signature(f="function"),
          definition=function (f, proto = NULL, ...) {
            if (!is.null(proto)) {
              if (!is.call(proto))
                stop(sQuote("proto")," must be an unevaluated call")
              prototype <- as.character(proto)
              fname <- prototype[1]
              args <- prototype[-1]
              if (is.function(f)&&(!all(args%in%names(formals(f)))))
                stop(sQuote(fname)," must be a function of prototype ",
                     deparse(proto),call.=FALSE)
            }
            new("pomp.fun",R.fun=f,mode=pompfunmode$Rfun)
          }
          )

setMethod(
          "pomp.fun",
          signature=signature(f="character"),
          definition=function (f, PACKAGE = NULL,
            obsnames = character(0), statenames = character(0),
            paramnames = character(0), covarnames = character(0), ...) {
            new(
                "pomp.fun",
                native.fun=f,
                PACKAGE=as.character(PACKAGE),
                mode=pompfunmode$native,
                obsnames=obsnames,
                statenames=statenames,
                paramnames=paramnames,
                covarnames=covarnames
                )
          }
          )

setMethod(
          "pomp.fun",
          signature=signature(f="Csnippet"),
          definition=function (f, PACKAGE, slotname = NULL, libname = NULL, 
            obsnames = character(0), statenames = character(0),
            paramnames = character(0), covarnames = character(0), ...) {
            if (is.null(slotname))
              stop("pomp.fun error: unspecified",sQuote("slotname"))
            if (is.null(libname))
              stop("pomp.fun error: unspecified",sQuote("libname"))
            slotname <- as.character(slotname)
            libname <- as.character(libname)            
            new(
                "pomp.fun",
                native.fun=render(fnames[[slotname]],name=libname),
                PACKAGE=as.character(libname),
                mode=pompfunmode$regNative,
                obsnames=obsnames,
                statenames=statenames,
                paramnames=paramnames,
                covarnames=covarnames
                )
          }
          )

setMethod(
          "pomp.fun",
          signature=signature(f="pomp.fun"),
          definition=function (f, ...) f
          )

setMethod(
          'show',
          signature=signature('pomp.fun'),
          definition=function (object) {
            mode <- object@mode
            if (mode==pompfunmode$Rfun) { # R function
              show(object@R.fun)
            } else if (mode==pompfunmode$native) { # user supplied native code
              cat("native function ",sQuote(object@native.fun),sep="")
              if (length(object@PACKAGE)>0)
                cat(", dynamically loaded from ",sQuote(object@PACKAGE),sep="")
              cat ("\n")
            } else if (mode==pompfunmode$regNative) { # built from Csnippets
              cat("native function ",sQuote(object@native.fun),sep="")
              if (length(object@PACKAGE)>0)
                cat(", dynamically loaded from ",sQuote(object@PACKAGE),sep="")
              cat ("\n")
            } else {
              cat("function not specified\n")
            }
          }
          )

setMethod(
          'print',
          'pomp.fun',
          function (x, ...) show(x)
          )
