## a class for functions that may be defined in R or using native routines
setClass(
         'pomp.fun',
         representation(
                        R.fun = 'function',
                        native.fun = 'character',
                        PACKAGE = 'character',
                        mode = 'integer',
                        address = 'externalptr'
                        ),
         prototype=prototype(
           R.fun=function(...)stop("unreachable error: please report this bug!"),
           native.fun=character(0),
           PACKAGE="",
           mode=-1L ## undefined
           )
         )

## constructor
pomp.fun <- function (f = NULL, PACKAGE, proto = NULL) {
  if (missing(PACKAGE)) PACKAGE <- ""
  if (!is.null(proto)) {
    if (!is.call(proto))
      stop(sQuote("proto")," must be an unevaluated call")
    prototype <- as.character(proto)
    fname <- prototype[1]
    args <- prototype[-1]
    if (is.function(f)&&(!all(args%in%names(formals(f)))))
      stop(sQuote(fname)," must be a function of prototype ",deparse(proto),call.=FALSE)
  }
  if (is(f,"pomp.fun")) {
    retval <- f
  } else if (is.function(f)) {
    retval <- new(
                  "pomp.fun",
                  R.fun=f,
                  mode=1L
                  )
  } else if (is.character(f)) {
    retval <- new(
                  "pomp.fun",
                  native.fun=f,
                  PACKAGE=PACKAGE,
                  mode=2L
                  )
  } else {
    retval <- new("pomp.fun")
  }
  retval
}

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
