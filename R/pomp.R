## This file defines 'pomp', the basic constructor of the pomp class

pomp <- function (data, times, t0, ..., rinit, rprocess, dprocess,
  rmeasure, dmeasure, skeleton, rprior, dprior, partrans,
  params, covar, zeronames,
  obsnames, statenames, paramnames, covarnames,
  PACKAGE, globals, cdir, cfile, shlib.args) {

  ep <- paste0("in ",sQuote("pomp"),": ")

  if (missing(data))
    stop(ep,sQuote("data")," is a required argument.",call.=FALSE)

  if (!inherits(data,what=c("data.frame","pomp","NULL")))
    stop(ep,sQuote("data")," must be a data frame or an object of class ",
      sQuote("pomp"),".",call.=FALSE)

  ## return as quickly as possible if no work is to be done
  if (nargs()==1) return(data)

  construct_pomp(
    data=data,times=times,t0=t0,...,
    rinit=rinit,rprocess=rprocess,dprocess=dprocess,
    rmeasure=rmeasure,dmeasure=dmeasure,
    skeleton=skeleton,rprior=rprior,dprior=dprior,partrans=partrans,
    params=params, covar=covar,zeronames=zeronames,
    obsnames=obsnames,statenames=statenames,paramnames=paramnames,
    covarnames=covarnames,PACKAGE=PACKAGE,
    globals=globals,cdir=cdir,cfile=cfile,shlib.args=shlib.args
  )
}

setMethod(
  "construct_pomp",
  signature=signature(data="data.frame"),
  definition = function (data, times, ...) {

    ep <- paste0("in ",sQuote("pomp"),": ")  # error prefix

    if (anyDuplicated(names(data))) {
      stop(ep,"names of data variables must be unique.", call.=FALSE)
    }
    if (missing(times))
      stop(ep,sQuote("times")," is a required argument",call.=FALSE)
    if ((is.numeric(times) && (times<1 || times>ncol(data) ||
        times!=as.integer(times))) ||
        (is.character(times) && (!(times%in%names(data)))) ||
        (!is.numeric(times) && !is.character(times)) || length(times)!=1) {
      stop(ep,"when ",sQuote("data")," is a data frame, ",sQuote("times"),
        " must identify a single column of ",sQuote("data"),
        " either by name or by index.",call.=FALSE)
    }
    if (is.numeric(times)) {
      tpos <- as.integer(times)
      timename <- names(data)[tpos]
    } else if (is.character(times)) {
      tpos <- match(times,names(data))
      timename <- times
    }
    times <- data[[tpos]]
    data <- do.call(rbind,lapply(data[-tpos],as.double))

    construct_pomp(data=data,times=times,...,timename=timename)
  }
)

setMethod(
  "construct_pomp",
  signature=signature(data="NULL"),
  definition = function (data, times, ...) {

    ep <- paste0("in ",sQuote("pomp"),": ")

    if (missing(times))
      stop(ep,sQuote("times")," is a required argument.",call.=FALSE)

    construct_pomp(data=array(dim=c(0,length(times))),times=times,...)
  }
)

setMethod(
  "construct_pomp",
  signature=signature(data="array"),
  definition = function (data, times, t0, ...,
    rinit, rprocess, dprocess, rmeasure, dmeasure, skeleton, rprior, dprior,
    partrans, params, covar, timename) {

    ep <- paste0("in ",sQuote("pomp"),": ")

    if (missing(rinit)) rinit <- NULL

    if (missing(rprocess) || is.null(rprocess)) {
      rprocess <- rproc_plugin()
    }

    if (missing(dprocess)) dprocess <- NULL
    if (missing(rmeasure)) rmeasure <- NULL
    if (missing(dmeasure)) dmeasure <- NULL

    if (missing(skeleton) || is.null(skeleton)) {
      skeleton <- skel_plugin()
    }

    if (missing(partrans) || is.null(partrans)) {
      partrans <- parameter_trans()
    }

    if (missing(rprior)) rprior <- NULL
    if (missing(dprior)) dprior <- NULL

    if (missing(params)) params <- numeric(0)
    if (missing(covar)) covar <- covariate_table()

    tryCatch(
      pomp.internal(
        data=data,
        times=times,
        t0=t0,
        timename=timename,
        rinit=rinit,
        rprocess=rprocess,
        dprocess=dprocess,
        rmeasure=rmeasure,
        dmeasure=dmeasure,
        skeleton=skeleton,
        dprior=dprior,
        rprior=rprior,
        partrans=partrans,
        params=params,
        covar=covar,
        ...
      ),
      error = function (e) {
        stop(ep,conditionMessage(e),call.=FALSE)
      }
    )
  }
)

setMethod(
  "construct_pomp",
  signature=signature(data="pomp"),
  definition = function (data, times, t0, ...,
    rinit, rprocess, dprocess, rmeasure, dmeasure, skeleton, rprior, dprior,
    partrans, params, covar, zeronames, timename) {

    ep <- paste0("in ",sQuote("pomp"),": ")  # error prefix

    if (missing(times)) {
      times <- data@times
    } else {
      time(data) <- times
    }

    if (missing(t0)) t0 <- data@t0
    if (missing(timename)) timename <- data@timename

    if (missing(rinit)) rinit <- data@rinit

    if (missing(rprocess)) {
      rprocess <- data@rprocess
    } else if (is.null(rprocess)) {
      rprocess <- rproc_plugin()
    }

    if (missing(dprocess)) dprocess <- data@dprocess
    if (missing(rmeasure)) rmeasure <- data@rmeasure
    if (missing(dmeasure)) dmeasure <- data@dmeasure

    if (missing(skeleton)) {
      skeleton <- data@skeleton
    } else if (is.null(skeleton)) {
      skeleton <- skel_plugin()
    }

    if (missing(partrans)) {
      partrans <- data@partrans
    } else if (is.null(partrans)) {
      partrans <- parameter_trans()
    }

    if (missing(rprior)) rprior <- data@rprior
    if (missing(dprior)) dprior <- data@dprior

    if (missing(params)) params <- data@params
    if (missing(covar)) covar <- data@covar
    if (missing(zeronames)) zeronames <- data@zeronames

    tryCatch(
      pomp.internal(
        data=data@data,
        times=times,
        t0=t0,
        timename=timename,
        rinit=rinit,
        rprocess=rprocess,
        dprocess=dprocess,
        rmeasure=rmeasure,
        dmeasure=dmeasure,
        skeleton=skeleton,
        rprior=rprior,
        dprior=dprior,
        partrans=partrans,
        covar=covar,
        zeronames=zeronames,
        params=params,
        .solibs=data@solibs,
        .userdata=data@userdata,
        ...
      ),
      error = function (e) {
        stop(ep,conditionMessage(e),call.=FALSE)
      }
    )
  }
)

pomp.internal <- function (data, times, t0, timename, ...,
  rinit, rprocess, dprocess, rmeasure, dmeasure, skeleton, rprior, dprior,
  partrans, params, covar, zeronames, obsnames, statenames,
  paramnames, covarnames, PACKAGE, globals, cdir, cfile, shlib.args,
  .userdata, .solibs = list(), verbose = getOption("verbose",FALSE)) {

  ep <- character(0)
  wp <- paste0("in ",sQuote("pomp"),": ")

  ## check times
  if (missing(times) || !is.numeric(times) || !all(is.finite(times)) ||
      (length(times)>1 && !all(diff(times)>0)))
    stop(ep,sQuote("times")," must be specified as an increasing sequence ",
      "of numbers.",call.=FALSE)
  storage.mode(times) <- "double"

  ## check t0
  if (missing(t0) || !is.numeric(t0) || !is.finite(t0) ||
      length(t0) > 1 ||  t0 > times[1])
    stop(ep,sQuote("t0")," must be specified as a single number not ",
      "greater than ",sQuote("times[1]"),".",call.=FALSE)
  storage.mode(t0) <- "double"

  if (missing(timename) || is.null(timename))
    timename <- "time"
  else
    timename <- as.character(timename)

  if (missing(.userdata)) .userdata <- list()
  added.userdata <- list(...)
  if (length(added.userdata)>0) {
    message(wp,"the following unrecognized argument(s) ",
      "will be stored for use by user-defined functions: ",
      paste(sQuote(names(added.userdata)),collapse=","))
    .userdata[names(added.userdata)] <- added.userdata
  }

  if (!is(rprocess,"rprocPlugin")) {
    stop(ep,sQuote("rprocess"),
      " must be specified using one of the plugins:\n",
      sQuote("onestep.sim"),", ",sQuote("discrete.time.sim"),
      ", ",sQuote("euler.sim"),", ",sQuote("gillespie.sim"),
      ", or ",sQuote("gillespie.hl.sim"),".",call.=FALSE)
  }

  if (!is(skeleton,"skelPlugin")) {
    stop(ep,sQuote("skeleton")," must be specified using either ",
      sQuote("map")," or ",sQuote("vectorfield"),".",call.=FALSE)
  }

  if (!is(partrans,"partransPlugin")) {
    stop(ep,sQuote("partrans")," must be specified using ",
      sQuote("parameter_trans"),".",call.=FALSE)
  }

  if (missing(statenames)) statenames <- NULL
  if (missing(paramnames)) paramnames <- NULL
  if (missing(obsnames)) obsnames <- NULL
  if (missing(covarnames)) covarnames <- NULL

  statenames <- as.character(statenames)
  paramnames <- as.character(paramnames)
  obsnames <- as.character(obsnames)
  covarnames <- as.character(covarnames)

  if (missing(zeronames)) zeronames <- NULL
  zeronames <- unique(as.character(zeronames))

  ## store the data as double-precision matrix
  storage.mode(data) <- "double"
  if (length(obsnames) == 0) obsnames <- rownames(data)

  ## check the parameters and force them to be double-precision
  if (length(params)>0) {
    if (is.null(names(params)) || !is.numeric(params) || any(names(params)==""))
      stop(sQuote("params")," must be a named numeric vector.",call.=FALSE)
  }
  storage.mode(params) <- "double"


  ## use default rinit?
  default.init <- is.null(rinit) ||
    (is(rinit,"pomp_fun") && rinit@mode == pompfunmode$undef )
  if (default.init) rinit <- pomp_fun(slotname="rinit")

  if (is(rinit,"Csnippet") && length(statenames)==0) {
    stop(ep,"when ",sQuote("rinit")," is provided as a C snippet, ",
      "you must also provide ",sQuote("statenames"),".",call.=FALSE)
  }

  ## by default, use flat improper prior
  if (is.null(dprior))
    dprior <- pomp_fun(f="_pomp_default_dprior",PACKAGE="pomp")

  ## check and arrange covariates
  if (is.null(covar)) {
    covar <- covariate_table()
  } else if (!is(covar,"covartable")) {
    stop(ep,"bad option for ",sQuote("covar"),".",call.=FALSE)
  }

  if (length(covarnames) == 0) {
    covarnames <- get_covariate_names(covar)
  } else {
    covar <- tryCatch(
      select_covariates(covar,covarnames),
      error = function (e) {
        stop(ep,conditionMessage(e),call.=FALSE)
      }
    )
  }

  hitches <- hitch(
    rinit=rinit,
    step.fn=rprocess@step.fn,
    rate.fn=rprocess@rate.fn,
    dprocess=dprocess,
    rmeasure=rmeasure,
    dmeasure=dmeasure,
    rprior=rprior,
    dprior=dprior,
    toEst=partrans@to,
    fromEst=partrans@from,
    skeleton=skeleton@skel.fn,
    templates=workhorse_templates,
    obsnames=obsnames,
    statenames=statenames,
    paramnames=paramnames,
    covarnames=covarnames,
    PACKAGE=PACKAGE,
    globals=globals,
    cfile=cfile,
    cdir=cdir,
    shlib.args=shlib.args,
    verbose=verbose
  )

  ## check to make sure 'covars' is included as an argument where needed
  covar_fun_warning(covar,"rprocess",hitches$funs,wp)
  covar_fun_warning(covar,"dprocess",hitches$funs,wp)
  covar_fun_warning(covar,"rmeasure",hitches$funs,wp)
  covar_fun_warning(covar,"dmeasure",hitches$funs,wp)
  covar_fun_warning(covar,"skeleton",hitches$funs,wp)

  ## check to see if covariate times embrace the data times
  covar_time_warning(covar,times,t0,wp)

  new(
    "pomp",
    data = data,
    times = times,
    t0 = t0,
    timename = timename,
    default.init = default.init,
    rinit = hitches$funs$rinit,
    rprocess = rproc_plugin(
      rprocess,
      step.fn=hitches$funs$step.fn,
      rate.fn=hitches$funs$rate.fn
    ),
    dprocess = hitches$funs$dprocess,
    dmeasure = hitches$funs$dmeasure,
    rmeasure = hitches$funs$rmeasure,
    skeleton = skel_plugin(
      skeleton,
      skel.fn=hitches$funs$skeleton
    ),
    dprior = hitches$funs$dprior,
    rprior = hitches$funs$rprior,
    partrans = parameter_trans(
      toEst=hitches$funs$toEst,
      fromEst=hitches$funs$fromEst
    ),
    params = params,
    covar = covar,
    zeronames = zeronames,
    solibs = if (is.null(hitches$lib)) {
      .solibs
    } else {
      c(list(hitches$lib),.solibs)
    },
    userdata = .userdata
  )
}
