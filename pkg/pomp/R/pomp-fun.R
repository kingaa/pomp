## a class for functions that may be defined in R or using native routines
setClass(
         'pomp.fun',
         representation(
                        R.fun = 'function',
                        native.fun = 'character',
                        PACKAGE = 'character',
                        mode = 'integer'
                        )
         )

pomp_native_code <- 2L
pomp_R_function <- 1L
pomp_undef_mode <- -1L

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
                  native.fun=character(0),
                  PACKAGE=PACKAGE,
                  mode=pomp_R_function
                  )
  } else if (is.character(f)) {
    retval <- new(
                  "pomp.fun",
                  R.fun=function(...)stop("unreachable error: please report this bug!"),
                  native.fun=f,
                  PACKAGE=PACKAGE,
                  mode=pomp_native_code
                  )
  } else {
    retval <- new(
                  "pomp.fun",
                  R.fun=function(...)stop(sQuote(fname)," not specified"),
                  native.fun=character(0),
                  PACKAGE=PACKAGE,
                  mode=pomp_undef_mode
                  )
  }
  retval
}

setMethod(
          'show',
          'pomp.fun',
          function (object) {
            mode <- object@mode
            if (mode==pomp_R_function) {
              show(object@R.fun)
            } else if (mode==pomp_native_code) {
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

get.pomp.fun <- function (object)
  .Call(get_pomp_fun,object)


## get.pomp.fun <- function (object) {
##   mode <- object@mode
##   if (mode==1L) {
##     f <- object@R.fun
##   } else if (mode==2L) {
##     f <- getNativeSymbolInfo(name=object@native.fun,PACKAGE=object@PACKAGE)$address
##   } else {
##     stop("function not specified")
##   }
##   list(f,as.integer(mode-1))
## }
