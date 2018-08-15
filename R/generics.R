## constituent components of a 'pomp' object
setGeneric("dmeasure",function(object,...)standardGeneric("dmeasure"))
setGeneric("rmeasure",function(object,...)standardGeneric("rmeasure"))
setGeneric("dprocess",function(object,...)standardGeneric("dprocess"))
setGeneric("rprocess",function(object,...)standardGeneric("rprocess"))
setGeneric("dprior",function(object,...)standardGeneric("dprior"))
setGeneric("rprior",function(object,...)standardGeneric("rprior"))
setGeneric("rinit",function(object,...)standardGeneric("rinit"))
setGeneric("skeleton",function(object,...)standardGeneric("skeleton"))

## functions to extract or call the components of a "pomp" object
setGeneric("states",function(object,...)standardGeneric("states"))
setGeneric("partrans",function(object,...)standardGeneric("partrans"))
setGeneric("logEvidence",function(object,...)standardGeneric("logEvidence"))

## internals
setGeneric("pomp_fun",function(f,...)standardGeneric("pomp_fun"))
## dynamic loading and unloading
setGeneric("pompLoad",function(object,...)standardGeneric("pompLoad"))
setGeneric("pompUnload",function(object,...)standardGeneric("pompUnload"))
setGeneric("solibs<-",function(object,...,value)standardGeneric("solibs<-"))

setGeneric("cond.logLik",function(object,...)standardGeneric("cond.logLik"))
setGeneric("cond.logEvidence",
  function(object,...)standardGeneric("cond.logEvidence"))

## values of probes
setGeneric("probevals",
  function(object,...)standardGeneric("probevals"))
## samples from prior and posterior distributions
setGeneric("prior_samples",
  function(object,...)standardGeneric("prior_samples"))
setGeneric("posterior_samples",
  function(object,...)standardGeneric("posterior_samples"))

## helper for ABC and PMCMC
setGeneric("covmat",function(object,...)standardGeneric("covmat"))
