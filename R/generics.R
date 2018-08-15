## basic methods
setGeneric("print",function(x,...)standardGeneric("print"))
setGeneric("plot",function(x,y,...)standardGeneric("plot"))
setGeneric("summary",function(object,...)standardGeneric("summary"))
setGeneric("window",function(x,...)standardGeneric("window"))
setGeneric("spy",function(object,...)standardGeneric("spy"))
setGeneric("concat",function(...)standardGeneric("concat"))
## continue an iterative calculation
setGeneric("continue",function(object,...)standardGeneric("continue"))

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
setGeneric("obs",function(object,...)standardGeneric("obs"))
setGeneric("time",function(x,...)standardGeneric("time"))
setGeneric("time<-",function(object,...,value)standardGeneric("time<-"))
setGeneric("coef",function(object,...)standardGeneric("coef"))
setGeneric("coef<-",function(object,...,value)standardGeneric("coef<-"))
setGeneric("states",function(object,...)standardGeneric("states"))
setGeneric("timezero",function(object,...)standardGeneric("timezero"))
setGeneric("timezero<-",function(object,...,value)standardGeneric("timezero<-"))
setGeneric("partrans",function(object,...)standardGeneric("partrans"))
setGeneric("logLik",function(object,...)standardGeneric("logLik"))
setGeneric("logEvidence",function(object,...)standardGeneric("logEvidence"))
## history of an iterative calculation
setGeneric("traces",function(object,...)standardGeneric("traces"))
setGeneric("parameter_trans",
  function(toEst,fromEst,...)standardGeneric("parameter_trans"))
setGeneric("covariate_table",
  function(...,times)standardGeneric("covariate_table"))

## internals
setGeneric("pomp_fun",function(f,...)standardGeneric("pomp_fun"))
## dynamic loading and unloading
setGeneric("pompLoad",function(object,...)standardGeneric("pompLoad"))
setGeneric("pompUnload",function(object,...)standardGeneric("pompUnload"))
setGeneric("solibs<-",function(object,...,value)standardGeneric("solibs<-"))

setGeneric("pred.mean",function(object,...)standardGeneric("pred.mean"))
setGeneric("pred.var",function(object,...)standardGeneric("pred.var"))
setGeneric("filter.mean",function(object,...)standardGeneric("filter.mean"))
setGeneric("filter.traj",function(object,...)standardGeneric("filter.traj"))
setGeneric("forecast",function(object,...)standardGeneric("forecast"))
setGeneric("cond.logLik",function(object,...)standardGeneric("cond.logLik"))
setGeneric("cond.logEvidence",
  function(object,...)standardGeneric("cond.logEvidence"))
setGeneric("eff.sample.size",
  function(object,...)standardGeneric("eff.sample.size"))

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

setMethod(
  "continue",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("continue"),": ",sQuote("object"),
      " is a required argument",call.=FALSE)
  }
)


setMethod(
  "continue",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("continue")," is not defined for objects of class ",
      sQuote(class(object)),call.=FALSE)
  }
)

setMethod(
  "filter.mean",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("filter.mean")," is not defined for objects of class ",
      sQuote(class(object)),call.=FALSE)
  }
)

setMethod(
  "filter.traj",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("filter.traj")," is not defined for objects of class ",
      sQuote(class(object)),call.=FALSE)
  }
)

setMethod(
  "pred.mean",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("pred.mean")," is not defined for objects of class ",
      sQuote(class(object)),call.=FALSE)
  }
)

setMethod(
  "pred.var",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("pred.var")," is not defined for objects of class ",
      sQuote(class(object)),call.=FALSE)
  }
)

setMethod(
  "eff.sample.size",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("eff.sample.size")," is not defined for objects of class ",
      sQuote(class(object)),call.=FALSE)
  }
)

setMethod(
  "traces",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("traces")," is not defined for objects of class ",
      sQuote(class(object)),call.=FALSE)
  }
)

setMethod(
  "concat",
  signature=signature(...="ANY"),
  definition=function(...) {
    stop(sQuote("c")," is not defined for objects of mixed class.",
      call.=FALSE)
  }
)
