##' The basic pomp class
##'
##' The basic class implementing a \acronym{POMP} model with data
##'
##' @rdname pomp_class
##' @include pomp_fun.R
##' @include csnippet.R safecall.R
##' @include rprocess_spec.R skeleton_spec.R
##' @include covariate_table.R parameter_trans.R
##' @keywords internal

setClass(
  "pomp",
  slots=c(
    data = "array",
    times = "numeric",
    t0 = "numeric",
    timename = "character",
    default.init = "logical",
    rinit = "pomp_fun",
    rprocess = "rprocPlugin",
    dprocess = "pomp_fun",
    dmeasure = "pomp_fun",
    rmeasure = "pomp_fun",
    dprior = "pomp_fun",
    rprior = "pomp_fun",
    skeleton = "skelPlugin",
    partrans = "partransPlugin",
    states = "array",
    params = "numeric",
    covar = "covartable",
    zeronames = "character",
    solibs = "list",
    userdata = "list"
  ),
  prototype=prototype(
    data=array(data=numeric(0),dim=c(0,0)),
    times=numeric(0),
    t0=numeric(0),
    timename="time",
    default.init=TRUE,
    rinit=pomp_fun(slotname="rinit"),
    rprocess=rproc_plugin(),
    dprocess=pomp_fun(slotname="dprocess"),
    dmeasure=pomp_fun(slotname="dmeasure"),
    rmeasure=pomp_fun(slotname="rmeasure"),
    dprior=pomp_fun(slotname="dprior"),
    rprior=pomp_fun(slotname="rprior"),
    skeleton=skel_plugin(),
    partrans=parameter_trans(),
    states=array(data=numeric(0),dim=c(0,0)),
    params=numeric(0),
    covar=covariate_table(),
    zeronames=character(0),
    solibs=list(),
    userdata=list()
  ),
  validity=function (object) {
    retval <- character(0)
    if (length(object@times)<1)
      retval <- append(retval,paste(sQuote("times"),"is a required argument"))
    if (!is.numeric(object@params) || (length(object@params)>0 && is.null(names(object@params))))
      retval <- append(retval,paste(sQuote("params"),"must be a named numeric vector"))
    if (ncol(object@data)!=length(object@times))
      retval <- append(retval,paste("the length of",sQuote("times"),"should match the number of observations"))
    if (length(object@t0)<1)
      retval <- append(retval,paste(sQuote("t0"),"is a required argument"))
    if (!is.numeric(object@t0) || !is.finite(object@t0) || length(object@t0)>1)
      retval <- append(retval,paste(sQuote("t0"),"must be a single number"))
    if (object@t0 > object@times[1])
      retval <- append(retval,paste("the zero-time",sQuote("t0"),
        "must occur no later than the first observation"))
    if (length(retval)==0) TRUE else retval
  }
)
