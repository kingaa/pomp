## a class for functions that may be defined in R or using native routines
setClass(
         'pomp.fun',
         representation(
                        R.fun = 'function',
                        native.fun = 'character',
                        PACKAGE = 'character',
                        use = 'integer'
                        )
         )

## constructor
pomp.fun <- function (f = NULL, PACKAGE, proto = NULL) {
  if (missing(PACKAGE))
    PACKAGE <- character(0)
  if (is.function(f)) {
    if (!is.null(proto)) {
      prototype <- strsplit(as.character(proto),split="\\(|\\)|\\,")
      fname <- prototype[1]
      args <- prototype[-1]
      if (!all(args%in%names(formals(f))))
        stop(sQuote(fname)," must be a function of prototype ",deparse(proto),call.=FALSE)
    }
    retval <- new(
                  "pomp.fun",
                  R.fun=f,
                  native.fun=character(0),
                  PACKAGE=PACKAGE,
                  use=1L
                  )
  } else if (is.character(f)) {
    retval <- new(
                  "pomp.fun",
                  R.fun=function(...)stop("unreachable error: please report this bug!"),
                  native.fun=f,
                  PACKAGE=PACKAGE,
                  use=2L
                  )
  } else {
    retval <- new(
                  "pomp.fun",
                  R.fun=function(...)stop(sQuote(fname)," not specified"),
                  native.fun=character(0),
                  PACKAGE=PACKAGE,
                  use=1L
                  )
  }
  retval
}

setMethod(
          'show',
          'pomp.fun',
          function (object) {
            if (object@use==1) {
              show(object@R.fun)
            } else {
              cat("native function ",sQuote(object@native.fun),sep="")
              if (length(object@PACKAGE)>0)
                cat(", dynamically loaded from ",sQuote(object@PACKAGE),sep="")
              cat ("\n")
            }
          }
          )

setMethod(
          'print',
          'pomp.fun',
          function (x, ...) show(x)
          )
