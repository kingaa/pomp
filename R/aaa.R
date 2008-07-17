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
                        skeleton = 'pomp.fun',
                        initializer = 'function',
                        states = 'array',
                        params = 'numeric',
                        covar = 'matrix',
                        tcovar = 'numeric',
                        obsindex = 'integer',
                        statenames = 'character',
                        paramnames = 'character',
                        covarnames = 'character',
                        PACKAGE = 'character',
                        userdata = 'list'
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

skeleton <- function (object, x, t, params, ...)
  stop("function 'skeleton' is undefined for objects of class '",class(object),"'")
setGeneric('skeleton')  

init.state <- function (object, params, t0, ...)
  stop("function 'init.state' is undefined for objects of class '",class(object),"'")
setGeneric('init.state')  

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

'coef<-' <- function (object, pars, ..., value)
  stop("function 'coef<' is undefined for objects of class '",class(object),"'")
setGeneric('coef<-')

states <- function (object, ...)
  stop("function 'states' is undefined for objects of class '",class(object),"'")
setGeneric('states')

bsmc <- function (object, ...)
  stop("function 'bsmc' is undefined for objects of class '",class(object),"'")
setGeneric('bsmc')

#lpspomp and pmpomp class
lpspomp <- function (object, ...)
  stop("function 'lpspomp' is undefined for objects of class '",class(object),"'")
setGeneric('lpspomp')

pmpomp <- function (object, ...)
  stop("function 'pmpomp' is undefined for objects of class '",class(object),"'")
setGeneric('pmpomp')

distance.hist <- function (object, ...)
  stop("function 'distance.hist' is undefined for objects of class '",class(object),"'")
setGeneric('distance.hist')

r.hist <- function (object, ...)
  stop("function 'r.hist' is undefined for objects of class '",class(object),"'")
setGeneric('r.hist')

sensitivity.spectrum.fit <- function (object, ...)
  stop("function 'sensitivity.spectrum.fit' is undefined for objects of class '",class(object),"'")
setGeneric('sensitivity.spectrum.fit')

# The class "lpspomp" (log power spectrum pomp). Its parts are:
# 1) A pomp, with filled "params" and "initializer" slots
# 2) Log power spectra of simulated time series
# 3) Log power spectra of data time series
# 4) The p-value for the "spectrum distance fit test"
# 5) The p-value for the "spectrum shape fit test"
# 6) The kernel width used for spectrum smoothing
# 7) Specifies the scale to which simulations and data are converted before taking power spectra
# 8) Specifies the kind of detrending that will be done on all time series before taking spectra
setClass("lpspomp",
         representation("pomp",
                        sps="array",
                        dps="array",
                        L2.distance.fit="numeric",
                        r.fit="numeric",
                        kernel.width="numeric",
                        tsscale="character",
                        detrend="numeric"),
        )

# The class "pmpomp." Its parts are:
# 1) A pomp, with certain slots filled.
# 2) A named list of probes. These should be R functions.
# 3) A list of optional arguments to probes.
# 4) A matrix for which each column has probe values for simualtions.
# Columns correspond to the probes listed in 4, and rows correspond
# to simulations.
# 5) A vector with entries the probe values on the data.
# 6) A vector with entries, for each probe, the fraction of simulations
# for which the probe is less than the value of the probe on the data.
# 7) Same but for fractions greater than.
setClass("pmpomp",
         representation("pomp",
                        probes="list",
                        probe.opt.args="list",
                        simvals="array",
                        datvals="numeric",
                        fractions.less="numeric",
                        fractions.more="numeric"),
        )

