## as of version 0.37-1 'pomp' is a generic function
setGeneric("pomp",function(data,...)standardGeneric("pomp"))

## this is the initial-condition setting function that is used by default
## it simply finds all parameters in the vector 'params' that have a name ending in '.0'
## and returns a vector with their values with names stripped of '.0'
default.initializer <- function (params, t0, ...) {
  ivpnames <- grep("\\.0$",names(params),value=TRUE)
  if (length(ivpnames)<1)
    stop("default initializer error: no parameter names ending in ",
         sQuote(".0")," found: see ",sQuote("pomp")," documentation")
  x <- params[ivpnames]
  names(x) <- sub("\\.0$","",ivpnames)
  x
}

## define the pomp class
setClass(
         'pomp',
         representation=representation(
           data = 'array',
           times = 'numeric',
           t0 = 'numeric',
           rprocess = 'function',
           dprocess = 'function',
           dmeasure = 'pomp.fun',
           rmeasure = 'pomp.fun',
           skeleton.type = 'character',
           skeleton = 'pomp.fun',
           skelmap.delta.t = 'numeric',
           initializer = 'function',
           states = 'array',
           params = 'numeric',
           covar = 'matrix',
           tcovar = 'numeric',
           obsnames = 'character',
           statenames = 'character',
           paramnames = 'character',
           covarnames = 'character',
           zeronames = 'character',
           has.trans = 'logical',
           par.trans = 'pomp.fun',
           par.untrans = 'pomp.fun',
           PACKAGE = 'character',
           userdata = 'list'
           )
         )
