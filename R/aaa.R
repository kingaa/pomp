## all definitions of classes and generics should go in this file!
## the partially-observed Markov process ("pomp") object

## define the pomp class
setClass(
         'pomp',
         representation(
                        data = 'array',
                        times = 'numeric',
                        t0 = 'numeric',
                        rprocess = 'function',
                        dprocess = 'function',
                        dmeasure = 'function',
                        rmeasure = 'function',
                        userdata = 'list'
                        )
         )

## define the mif class
setClass(
         'mif',
         representation(
                        'pomp',
                        ivps = 'character',
                        pars = 'character',
                        Nmif = 'integer',
                        particles = 'function',
                        alg.pars = 'list',
                        coef = 'numeric',
                        random.walk.sd = 'numeric',
                        pred.mean = 'matrix',
                        pred.var = 'matrix',
                        filter.mean = 'matrix',
                        conv.rec = 'matrix',
                        eff.sample.size = 'numeric',
                        cond.loglik = 'numeric',
                        loglik = 'numeric'
                        )
         )

## functions to extract or call the components of a "pomp" object
data.array <- function (object, ...)
  stop("function 'data.array' is undefined for objects of class '",class(object),"'")
setGeneric('data.array')  

rprocess <- function (object, xstart, times, params, ...)
  stop("function 'rprocess' is undefined for objects of class '",class(object),"'")
setGeneric('rprocess')  

dprocess <- function (object, x, times, params, log = FALSE, ...)
  stop("function 'dprocess' is undefined for objects of class '",class(object),"'")
setGeneric('dprocess')  

rmeasure <- function (object, x, times, params, ...)
  stop("function 'rmeasure' is undefined for objects of class '",class(object),"'")
setGeneric('rmeasure')  

dmeasure <- function (object, y, x, times, params, log = FALSE, ...)
  stop("function 'dmeasure' is undefined for objects of class '",class(object),"'")
setGeneric('dmeasure')  

## particle filter
pfilter <- function (object, ...)
  stop("function 'pfilter' is undefined for objects of class '",class(object),"'")
setGeneric('pfilter')  

## MIF algorithm functions
mif <- function (object, ... )
  stop("function 'mif' is undefined for objects of class '",class(object),"'")
setGeneric('mif')

particles <- function (object, ...)
  stop("function 'particles' is undefined for objects of class '",class(object),"'")
setGeneric('particles')  

pred.mean <- function (object, ...)
  stop("function 'pred.mean' is undefined for objects of class '",class(object),"'")
setGeneric('pred.mean')  

pred.var <- function (object, ...)
  stop("function 'pred.var' is undefined for objects of class '",class(object),"'")
setGeneric('pred.var')  

filter.mean <- function (object, ...)
  stop("function 'filter.mean' is undefined for objects of class '",class(object),"'")
setGeneric('filter.mean')  

conv.rec <- function (object, ...)
  stop("function 'conv.rec' is undefined for objects of class '",class(object),"'")
setGeneric('conv.rec')  

continue <- function (object, ... )
  stop("function 'continue' is undefined for objects of class '",class(object),"'")
setGeneric('continue')

