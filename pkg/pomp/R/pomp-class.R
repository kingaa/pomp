## this is the initial-condition setting function that is used by default
## it simply finds all parameters in the vector 'params' that have a name ending in '.0'
## and returns a vector with their values with names stripped of '.0'
default.initializer <- function (params, t0, ...) {
  ivpnames <- grep("\\.0$",names(params),value=TRUE)
  if (length(ivpnames)<1)
    stop("default initializer error: no parameter names ending in ",
         sQuote(".0")," found: see ",sQuote("pomp")," documentation")
  setNames(params[ivpnames],sub("\\.0$","",ivpnames))
}

## define the pomp class
setClass(
         'pomp',
         slots=c(
           data = 'array',
           times = 'numeric',
           t0 = 'numeric',
           rprocess = 'function',
           dprocess = 'function',
           dmeasure = 'pomp.fun',
           rmeasure = 'pomp.fun',
           dprior = 'pomp.fun',
           rprior = 'pomp.fun',
           skeleton.type = 'character',
           skeleton = 'pomp.fun',
           skelmap.delta.t = 'numeric',
           initializer = 'function',
           states = 'array',
           params = 'numeric',
           covar = 'matrix',
           tcovar = 'numeric',
           zeronames = 'character',
           has.trans = 'logical',
           par.trans = 'pomp.fun',
           par.untrans = 'pomp.fun',
           solibfile = 'list',
           userdata = 'list'
           ),
         prototype=prototype(
           data=array(data=numeric(0),dim=c(0,0)),
           times=numeric(0),
           t0=numeric(0),
           rprocess=function(xstart,times,params,...)stop(sQuote("rprocess")," not specified"),
           dprocess=function(x,times,params,log=FALSE,...)stop(sQuote("dprocess")," not specified"),
           dmeasure=pomp.fun(),
           rmeasure=pomp.fun(),
           dprior=pomp.fun(),
           rprior=pomp.fun(),
           skeleton.type="map",
           skeleton=pomp.fun(),
           skelmap.delta.t=as.numeric(NA),
           initializer=default.initializer,
           states=array(data=numeric(0),dim=c(0,0)),
           params=numeric(0),
           covar=array(data=numeric(0),dim=c(0,0)),
           tcovar=numeric(0),
           zeronames=character(0),
           has.trans=FALSE,
           par.trans=pomp.fun(),
           par.untrans=pomp.fun(),
           solibfile=list(),
           userdata=list()
           ),
         validity=function (object) {
           retval <- character(0)
           if (length(object@data)<1)
             retval <- append(retval,paste(sQuote("data"),"is a required argument"))
           if (length(object@times)<1)
             retval <- append(retval,paste(sQuote("times"),"is a required argument"))
           if (!is.numeric(object@params) || (length(object@params)>0 && is.null(names(object@params))))
             retval <- append(retval,paste(sQuote("params"),"must be a named numeric vector"))
           if (ncol(object@data)!=length(object@times))
             retval <- append(retval,paste("the length of",sQuote("times"),"should match the number of observations"))
           if (length(object@t0)<1)
             retval <- append(retval,paste(sQuote("t0"),"is a required argument"))
           if (length(object@t0)>1)
             retval <- append(retval,paste(sQuote("t0"),"must be a single number"))
           if (object@t0 > object@times[1])
             retval <- append(retval,paste("the zero-time",sQuote("t0"),
                                           "must occur no later than the first observation"))
           if (object@skelmap.delta.t <= 0)
             retval <- append(retval,paste(sQuote("skelmap.delta.t"),"must be positive"))
           if (!all(c('xstart','times','params','...')%in%names(formals(object@rprocess))))
             retval <- append(
                              retval,
                              paste(
                                    sQuote("rprocess"),"must be a function of prototype",
                                    sQuote("rprocess(xstart,times,params,...)")
                                    )
                              )
           if (!all(c('x','times','params','log','...')%in%names(formals(object@dprocess))))
             retval <- append(
                              retval,
                              paste(
                                    sQuote("dprocess"),"must be a function of prototype",
                                    sQuote("dprocess(x,times,params,log,...)")
                                    )
                              )
           if (!all(c('params','t0','...')%in%names(formals(object@initializer))))
             retval <- append(
                              retval,
                              paste(
                                    sQuote("initializer"),"must be a function of prototype",
                                    sQuote("initializer(params,t0,...)")
                                    )
                              )
           if (length(object@tcovar)!=nrow(object@covar)) {
             retval <- append(
                              retval,
                              paste(
                                    "the length of",sQuote("tcovar"),
                                    "should match the number of rows of",sQuote("covar")
                                    )
                              )
           }
           if (!is.numeric(object@tcovar))
             retval <- append(
                              retval,
                              paste(
                                    sQuote("tcovar"),
                                    "must either be a numeric vector or must name a numeric vector in the data frame",
                                    sQuote("covar")
                                    )
                              )
           
           if (length(retval)==0) TRUE else retval
         }
         )
