## a class for functions that may be defined in R,
## using pre-written native routines, or C snippets

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
           mode=-1L, ## undefined behavior
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
            new("pomp.fun",R.fun=f,mode=1L)
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
                mode=2L,
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
                mode=2L,
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
          'pomp.fun',
          function (object) {
            mode <- object@mode
            if (mode==1L) {
              show(object@R.fun)
            } else if (mode==2L) {
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
