## all definitions of classes and generics should go in this file!
## the partially-observed Markov process ("pomp") object

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

## define the pomp class
setClass(
         'pomp',
         representation(
                        data = 'array',
                        times = 'numeric',
                        t0 = 'numeric',
                        rprocess = 'function',
                        dprocess = 'function',
                        dmeasure = 'pomp.fun',
                        rmeasure = 'pomp.fun',
                        skeleton.type = 'character',
                        skeleton = 'pomp.fun',
                        initializer = 'function',
                        states = 'array',
                        params = 'numeric',
                        covar = 'matrix',
                        tcovar = 'numeric',
                        obsnames = 'character',
                        statenames = 'character',
                        paramnames = 'character',
                        covarnames = 'character',
                        PACKAGE = 'character',
                        userdata = 'list',
                        call = "call"
                        )
         )

## define the mif class
setClass(
         'mif',
         representation(
                        ivps = 'character',
                        pars = 'character',
                        Nmif = 'integer',
                        particles = 'function',
                        alg.pars = 'list',
                        random.walk.sd = 'numeric',
                        pred.mean = 'matrix',
                        pred.var = 'matrix',
                        filter.mean = 'matrix',
                        conv.rec = 'matrix',
                        eff.sample.size = 'numeric',
                        cond.loglik = 'numeric',
                        loglik = 'numeric'
                        ),
         contains='pomp'
         )

## functions to extract or call the components of a "pomp" object
data.array <- function (object, ...)
  stop("function ",sQuote("data.array")," is undefined for objects of class ",sQuote(class(object)))
setGeneric('data.array')  

"time<-" <- function (object, ..., value)
  stop("function ",sQuote("time<-")," is undefined for objects of class ",sQuote(class(object)))
setGeneric("time<-")  

rprocess <- function (object, xstart, times, params, ...)
  stop("function ",sQuote("rprocess")," is undefined for objects of class ",sQuote(class(object)))
setGeneric('rprocess')  

dprocess <- function (object, x, times, params, log = FALSE, ...)
  stop("function ",sQuote("dprocess")," is undefined for objects of class ",sQuote(class(object)))
setGeneric('dprocess')  

rmeasure <- function (object, x, times, params, ...)
  stop("function ",sQuote("rmeasure")," is undefined for objects of class ",sQuote(class(object)))
setGeneric('rmeasure')  

dmeasure <- function (object, y, x, times, params, log = FALSE, ...)
  stop("function ",sQuote("dmeasure")," is undefined for objects of class ",sQuote(class(object)))
setGeneric('dmeasure')  

skeleton <- function (object, x, t, params, ...)
  stop("function ",sQuote("skeleton")," is undefined for objects of class ",sQuote(class(object)))
setGeneric('skeleton')  

init.state <- function (object, params, t0, ...)
  stop("function ",sQuote("init.state")," is undefined for objects of class ",sQuote(class(object)))
setGeneric('init.state')  

## particle filter
pfilter <- function (object, ...)
  stop("function ",sQuote("pfilter")," is undefined for objects of class ",sQuote(class(object)))
setGeneric('pfilter')  

## MIF algorithm functions
mif <- function (object, ... )
  stop("function ",sQuote("mif")," is undefined for objects of class ",sQuote(class(object)))
setGeneric('mif')

particles <- function (object, ...)
  stop("function ",sQuote("particles")," is undefined for objects of class ",sQuote(class(object)))
setGeneric('particles')  

pred.mean <- function (object, ...)
  stop("function ",sQuote("pred.mean")," is undefined for objects of class ",sQuote(class(object)))
setGeneric('pred.mean')  

pred.var <- function (object, ...)
  stop("function ",sQuote("pred.var")," is undefined for objects of class ",sQuote(class(object)))
setGeneric('pred.var')  

filter.mean <- function (object, ...)
  stop("function ",sQuote("filter.mean")," is undefined for objects of class ",sQuote(class(object)))
setGeneric('filter.mean')  

conv.rec <- function (object, ...)
  stop("function ",sQuote("conv.rec")," is undefined for objects of class ",sQuote(class(object)))
setGeneric('conv.rec')  

continue <- function (object, ... )
  stop("function ",sQuote("continue")," is undefined for objects of class ",sQuote(class(object)))
setGeneric('continue')

"coef<-" <- function (object, pars, ..., value)
  stop("function ",sQuote("coef<-")," is undefined for objects of class ",sQuote(class(object)))
setGeneric("coef<-")

states <- function (object, ...)
  stop("function ",sQuote("states")," is undefined for objects of class ",sQuote(class(object)))
setGeneric('states')

trajectory <- function (object, params, times, ...)
  stop("function ",sQuote("trajectory")," is undefined for objects of class ",sQuote(class(object)))
setGeneric('trajectory')
