## basic methods
setGeneric("construct_pomp",function(data,...)standardGeneric("construct_pomp"))
setGeneric("print",function(x,...)standardGeneric("print"))
setGeneric("plot",function(x,y,...)standardGeneric("plot"))
setGeneric("summary",function(object,...)standardGeneric("summary"))
setGeneric("window",function(x,...)standardGeneric("window"))
setGeneric("spy",function(object,...)standardGeneric("spy"))
setGeneric("concat",function(...)standardGeneric("concat"))
## continue an iterative calculation
setGeneric("continue",function(object,...)standardGeneric("continue"))
## history of an iterative calculation
setGeneric("traces",function(object,...)standardGeneric("traces"))


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
setGeneric("solibs<-",function(object,...,value)standardGeneric("solibs<-"))

## internals
setGeneric("pomp_fun",function(f,...)standardGeneric("pomp_fun"))

## prediction mean
setGeneric("pred.mean",function(object,...)standardGeneric("pred.mean"))
## prediction variance
setGeneric("pred.var",function(object,...)standardGeneric("pred.var"))
## filter mean
setGeneric("filter.mean",function(object,...)standardGeneric("filter.mean"))
## filter trajectory
setGeneric("filter.traj",function(object,...)standardGeneric("filter.traj"))
## conditional log likelihood
setGeneric("cond.logLik",function(object,...)standardGeneric("cond.logLik"))
setGeneric("cond.logEvidence",function(object,...)standardGeneric("cond.logEvidence"))
## effective sample size
setGeneric("eff.sample.size",function(object,...)standardGeneric("eff.sample.size"))
## values of probes
setGeneric("probevals",function(object,...)standardGeneric("probevals"))
## samples from prior and posterior distributions
setGeneric("prior_samples",function(object,...)standardGeneric("prior_samples"))
setGeneric("posterior_samples",function(object,...)standardGeneric("posterior_samples"))
## stochastic simulation
setGeneric("simulate",function(object,nsim=1,seed=NULL,...)standardGeneric("simulate"))

## deterministic trajectory computation
setGeneric("trajectory",function(object,...)standardGeneric("trajectory"))
## trajectory matching
setGeneric("traj.match.objfun",function(object,...)standardGeneric("traj.match.objfun"))
setGeneric("traj.match",function(object,...)standardGeneric("traj.match"))

## ABC algorithm functions
setGeneric('abc',function(object,...)standardGeneric("abc"))

## Kalman filter methods
setGeneric("enkf",function(object,...)standardGeneric("enkf"))
setGeneric("eakf",function(object,...)standardGeneric("eakf"))

## Bayesian SMC (Liu & West)
setGeneric("bsmc",function(object,...)standardGeneric("bsmc"))
setGeneric("bsmc2",function(object,...)standardGeneric("bsmc2"))

## basic SMC (particle filter)
setGeneric("pfilter",function(object,...)standardGeneric("pfilter"))

## particle Markov chain Monte Carlo (PMCMC)
setGeneric('pmcmc',function(object,...)standardGeneric("pmcmc"))

## helper for ABC and PMCMC
setGeneric("covmat",function(object,...)standardGeneric("covmat"))

## nonlinear forecasting
setGeneric('nlf',function(object,...)standardGeneric("nlf"))

## iterated filtering
setGeneric('mif',function(object,...)standardGeneric("mif"))
setGeneric("mif2",function(object,...)standardGeneric("mif2"))

## probe matching and synthetic likelihood
setGeneric("probe",function(object,probes,...)standardGeneric("probe"))
setGeneric("probe.match.objfun",function(object,...)standardGeneric("probe.match.objfun"))
setGeneric("probe.match",function(object,...)standardGeneric("probe.match"))

## power spectrum
setGeneric("spect",function(object,...)standardGeneric("spect"))
setGeneric("spect.match",function(object,...)standardGeneric("spect.match"))

## dynamic loading and unloading
setGeneric("pompLoad",function(object,...)standardGeneric("pompLoad"))
setGeneric("pompUnload",function(object,...)standardGeneric("pompUnload"))

setMethod(
  "traj.match",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("traj.match"),": ",sQuote("object")," is a required argument",call.=FALSE)
  }
)

setMethod(
  "traj.match",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("traj.match")," is not defined when ",sQuote("object")," is of class ",sQuote(class(object)),call.=FALSE)
  }
)

setMethod(
  "traj.match.objfun",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("traj.match.objfun"),": ",sQuote("object")," is a required argument",call.=FALSE)
  }
)

setMethod(
  "traj.match.objfun",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("traj.match.objfun")," is not defined when ",sQuote("object")," is of class ",sQuote(class(object)),call.=FALSE)
  }
)

setMethod(
  "abc",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("abc"),": ",sQuote("object")," is a required argument",call.=FALSE)
  }
)

setMethod(
  "abc",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("abc")," is not defined when ",sQuote("object")," is of class ",sQuote(class(object)),call.=FALSE)
  }
)

setMethod(
  "eakf",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("eakf"),": ",sQuote("object")," is a required argument",call.=FALSE)
  }
)

setMethod(
  "eakf",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("eakf")," is not defined when ",sQuote("object")," is of class ",sQuote(class(object)),call.=FALSE)
  }
)

setMethod(
  "enkf",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("enkf"),": ",sQuote("object")," is a required argument",call.=FALSE)
  }
)

setMethod(
  "enkf",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("enkf")," is not defined when ",sQuote("object")," is of class ",sQuote(class(object)),call.=FALSE)
  }
)

setMethod(
  "bsmc2",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("bsmc2"),": ",sQuote("object")," is a required argument",call.=FALSE)
  }
)

setMethod(
  "bsmc2",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("bsmc2")," is not defined when ",sQuote("object")," is of class ",sQuote(class(object)),call.=FALSE)
  }
)

setMethod(
  "pfilter",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("pfilter"),": ",sQuote("object")," is a required argument",call.=FALSE)
  }
)

setMethod(
  "pfilter",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("pfilter")," is not defined when ",sQuote("object")," is of class ",sQuote(class(object)),call.=FALSE)
  }
)

setMethod(
  "pmcmc",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("pmcmc"),": ",sQuote("object")," is a required argument",call.=FALSE)
  }
)

setMethod(
  "pmcmc",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("pmcmc")," is not defined when ",sQuote("object")," is of class ",sQuote(class(object)),call.=FALSE)
  }
)

setMethod(
  "nlf",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("nlf"),": ",sQuote("object")," is a required argument",call.=FALSE)
  }
)

setMethod(
  "nlf",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("nlf")," is not defined when ",sQuote("object")," is of class ",sQuote(class(object)),call.=FALSE)
  }
)

setMethod(
  "mif2",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("mif2"),": ",sQuote("object")," is a required argument",call.=FALSE)
  }
)

setMethod(
  "mif2",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("mif2")," is not defined when ",sQuote("object")," is of class ",sQuote(class(object)),call.=FALSE)
  }
)

setMethod(
  "probe",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("probe"),": ",sQuote("object")," is a required argument",call.=FALSE)
  }
)

setMethod(
  "probe",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("probe")," is not defined when ",sQuote("object")," is of class ",sQuote(class(object)),call.=FALSE)
  }
)

setMethod(
  "probe.match",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("probe.match"),": ",sQuote("object")," is a required argument",call.=FALSE)
  }
)

setMethod(
  "probe.match",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("probe.match")," is not defined when ",sQuote("object")," is of class ",sQuote(class(object)),call.=FALSE)
  }
)

setMethod(
  "probe.match.objfun",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("probe.match.objfun"),": ",sQuote("object")," is a required argument",call.=FALSE)
  }
)

setMethod(
  "probe.match.objfun",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("probe.match.objfun")," is not defined when ",sQuote("object")," is of class ",sQuote(class(object)),call.=FALSE)
  }
)

setMethod(
  "spect",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("spect"),": ",sQuote("object")," is a required argument",call.=FALSE)
  }
)

setMethod(
  "spect",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("spect")," is not defined when ",sQuote("object")," is of class ",sQuote(class(object)),call.=FALSE)
  }
)

setMethod(
  "continue",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("continue"),": ",sQuote("object")," is a required argument",call.=FALSE)
  }
)

setMethod(
  "continue",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("continue")," is not defined when ",sQuote("object")," is of class ",sQuote(class(object)),call.=FALSE)
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
