## This file defines 'pomp', the basic constructor of the pomp class

pomp <- function (data, times, t0, ..., rprocess, dprocess,
                  rmeasure, dmeasure, measurement.model,
                  skeleton, initializer, rprior, dprior, params, covar, tcovar,
                  obsnames, statenames, paramnames, covarnames, zeronames,
                  PACKAGE, fromEstimationScale, toEstimationScale,
                  globals, cdir, cfile, shlib.args) {

  ep <- paste0("in ",sQuote("pomp"),": ")  # error prefix

  if (missing(data))
    stop(ep,sQuote("data")," is a required argument.",call.=FALSE)

  ## return as quickly as possible if no work is to be done
  if (nargs()==1) return(data)

  if (!is(data,"data.frame") && !is(data,"pomp"))
    stop(ep,sQuote("data")," must be a data frame or an object of class ",
         sQuote("pomp"),call.=FALSE)

  ## if one transformation is supplied, then both must be
  fes <- missing(fromEstimationScale)
  tes <- missing(toEstimationScale)
  if (xor(fes,tes))
    stop(ep,"if one of ",sQuote("toEstimationScale"),", ",
         sQuote("fromEstimationScale")," is supplied, then so must the other",
         call.=FALSE)

  ## if 'measurement model' is specified as a formula, this overrides
  ## specification of 'rmeasure' or 'dmeasure'
  if (!missing(measurement.model)) {
    if (!(missing(dmeasure) || is.null(dmeasure)) ||
        !(missing(rmeasure) || is.null(rmeasure)))
      warning(ep,"specifying ",sQuote("measurement.model"),
              " overrides specification of ",
              sQuote("rmeasure")," and ",sQuote("dmeasure"),".",
              call.=FALSE)
    mm <- measform2pomp(measurement.model)
    rmeasure <- mm$rmeasure
    dmeasure <- mm$dmeasure
  }

  ## the deterministic skeleton involves 'skeleton', 'skeleton.type', and
  ## 'skelmap.delta.t'
  skel.type <- "undef"
  skelmap.delta.t <- 1
  if (missing(skeleton)) {
    skeleton <- NULL
  } else if (is.null(skeleton)) {
    skel.type <- "remove"
  } else if (is(skeleton,"safecall")) {
    skeleton <- skeleton@call
    flist <- list(
      map=function (f, delta.t = 1) {
        skel.type <<- "map"
        skelmap.delta.t <<- as.numeric(delta.t)
        if (skelmap.delta.t <= 0)
          stop("in ",sQuote("map"),", ",sQuote("delta.t"),
               " must be positive",call.=FALSE)
        f
      },
      vectorfield=function (f) {
        skel.type <<- "vectorfield"
        f
      }
    )
    skeleton <- eval(skeleton,envir=flist,enclos=parent.frame())
  } else {
    stop(ep,sQuote("skeleton")," must be specified as either a ",
         sQuote("vectorfield")," or a ",sQuote("map"),".",call.=FALSE)
  }

  construct_pomp(
    data=data,times=times,t0=t0,...,
    rprocess=rprocess,dprocess=dprocess,
    rmeasure=rmeasure,dmeasure=dmeasure,
    initializer=initializer,
    skel.type=skel.type,skelmap.delta.t=skelmap.delta.t,skeleton=skeleton,
    rprior=rprior,dprior=dprior,params=params,
    covar=covar,tcovar=tcovar,
    obsnames=obsnames,statenames=statenames,paramnames=paramnames,
    covarnames=covarnames,zeronames=zeronames,PACKAGE=PACKAGE,
    fromEstimationScale=fromEstimationScale,toEstimationScale=toEstimationScale,
    globals=globals,cdir=cdir,cfile=cfile,shlib.args=shlib.args
  )
}

setMethod("construct_pomp",
          signature=signature(data="data.frame"),
          definition = function (data, times, ...) {

            ep <- paste0("in ",sQuote("pomp"),": ")  # error prefix

            if (missing(times)) stop(ep,sQuote("times")," is a required argument",call.=FALSE)
            if ((is.numeric(times) && (times<1 || times>ncol(data) || times!=as.integer(times))) ||
                (is.character(times) && (!(times%in%names(data)))) ||
                (!is.numeric(times) && !is.character(times)) ||
                length(times)!=1) {
              stop(ep,"when ",sQuote("data")," is a data frame, ",sQuote("times"),
                   " must identify a single column of ",sQuote("data"),
                   " either by name or by index.",call.=FALSE)
            }
            if (is.numeric(times)) {
              tpos <- as.integer(times)
            } else if (is.character(times)) {
              tpos <- match(times,names(data))
            }
            times <- data[[tpos]]
            data <- do.call(rbind, lapply(data[-tpos], as.numeric))

            construct_pomp(data=data,times=times,...)
          }
)

setMethod("construct_pomp",
          signature=signature(data="missing"),
          definition = function (
            data, times, t0, ..., rprocess, rmeasure, obsnames
          ) {

            ep <- paste0("in ",sQuote("pomp"),": ")

            if (missing(times) || !is.numeric(times))
              stop(ep,sQuote("times")," must be supplied")
            if (missing(t0) || !is.numeric(t0))
              stop(ep,sQuote("t0")," must be supplied")
            if (missing(rprocess))
              stop(ep,sQuote("rprocess")," must be supplied")
            if (missing(rmeasure))
              stop(ep,sQuote("rmeasure")," must be supplied")
            if (missing(obsnames))
              stop(ep,sQuote("obsnames")," must be supplied")

            data <- array(dim=c(length(obsnames),length(times)),
                          dimnames=list(obsnames,NULL))

            construct_pomp(data=data,times=times,t0=t0,
                           rprocess=rprocess,rmeasure=rmeasure,obsnames=obsnames,...)
          }
)

setMethod("construct_pomp",
          signature=signature(data="array"),
          definition = function (data, times, t0, ...,
                                 rprocess, dprocess,
                                 rmeasure, dmeasure,
                                 skel.type, skelmap.delta.t, skeleton,
                                 initializer,
                                 rprior, dprior,
                                 params, covar, tcovar,
                                 fromEstimationScale, toEstimationScale) {

            ep <- paste0("in ",sQuote("pomp"),": ")  # error prefix

            if (missing(initializer)) initializer <- NULL
            if (missing(rprocess)) rprocess <- NULL
            if (missing(dprocess)) dprocess <- NULL
            if (missing(rmeasure)) rmeasure <- NULL
            if (missing(dmeasure)) dmeasure <- NULL
            if (missing(rprior)) rprior <- NULL
            if (missing(dprior)) dprior <- NULL
            if (missing(fromEstimationScale)) fromEstimationScale <- NULL
            if (missing(toEstimationScale)) toEstimationScale <- NULL

            if (missing(params)) params <- numeric(0)
            if (missing(covar)) covar <- NULL
            if (missing(tcovar)) tcovar <- NULL

            tryCatch(
              pomp.internal(
                data=data,
                times=times,
                t0=t0,
                rprocess=rprocess,
                dprocess=dprocess,
                rmeasure=rmeasure,
                dmeasure=dmeasure,
                dprior=dprior,
                rprior=rprior,
                skeleton=skeleton,
                skel.type=skel.type,
                skelmap.delta.t=skelmap.delta.t,
                initializer=initializer,
                params=params,
                covar=covar,
                tcovar=tcovar,
                fromEstimationScale=fromEstimationScale,
                toEstimationScale=toEstimationScale,
                ...
              ),
              error = function (e) {
                stop(ep,conditionMessage(e),call.=FALSE)
              }
            )
          }
)

setMethod("construct_pomp",
          signature=signature(data="pomp"),
          definition = function (
            data, times, t0, ...,
            rprocess, dprocess,
            rmeasure, dmeasure,
            skel.type, skelmap.delta.t, skeleton,
            initializer,
            rprior, dprior,
            params, covar, tcovar,
            zeronames,
            fromEstimationScale, toEstimationScale
          ) {

            ep <- paste0("in ",sQuote("pomp"),": ")  # error prefix

            if (missing(times)) {
              times <- data@times
            } else {
              time(data) <- times
            }

            if (missing(t0)) t0 <- data@t0

            if (missing(initializer)) initializer <- data@initializer
            if (missing(rprocess)) rprocess <- data@rprocess
            if (missing(dprocess)) dprocess <- data@dprocess
            if (missing(rmeasure)) rmeasure <- data@rmeasure
            if (missing(dmeasure)) dmeasure <- data@dmeasure
            if (missing(rprior)) rprior <- data@rprior
            if (missing(dprior)) dprior <- data@dprior
            if (missing(fromEstimationScale)) fromEstimationScale <- data@from.trans
            if (missing(toEstimationScale)) toEstimationScale <- data@to.trans

            if (missing(params)) params <- data@params
            if (missing(covar)) covar <- data@covar
            if (missing(tcovar)) tcovar <- data@tcovar
            if (missing(zeronames)) zeronames <- data@zeronames

            if (skel.type == "remove") {
              skel.type <- "undef"
              skelmap.delta.t <- 1
              skeleton <- NULL
            } else if (skel.type == "undef") {
              skel.type <- data@skeleton.type
              skelmap.delta.t <- data@skelmap.delta.t
              skeleton <- data@skeleton
            }

            tryCatch(
              pomp.internal(
                data=data@data,
                times=times,
                t0=t0,
                rprocess=rprocess,
                dprocess=dprocess,
                rmeasure=rmeasure,
                dmeasure=dmeasure,
                dprior=dprior,
                rprior=rprior,
                skeleton=skeleton,
                skel.type=skel.type,
                skelmap.delta.t=skelmap.delta.t,
                initializer=initializer,
                covar=covar,
                tcovar=tcovar,
                zeronames=zeronames,
                fromEstimationScale=fromEstimationScale,
                toEstimationScale=toEstimationScale,
                params=params,
                .solibs=data@solibs,
                userdata=data@userdata,
                ...
              ),
              error = function (e) {
                stop(ep,conditionMessage(e),call.=FALSE)
              }
            )
          }
)

setMethod("construct_pomp",
          signature=signature(data="ANY"),
          definition = function (data, ...) {
            ep <- paste0("in ",sQuote("pomp"),": ")
            stop(ep,sQuote("data")," must be a data frame or an object of class ",
                 sQuote("pomp"),call.=FALSE)
          }
)

pomp.internal <- function (
  data, times, t0, rprocess, dprocess,
  rmeasure, dmeasure,
  skeleton, skel.type, skelmap.delta.t,
  initializer, rprior, dprior,
  params, covar, tcovar,
  obsnames, statenames, paramnames, covarnames,
  zeronames, PACKAGE,
  fromEstimationScale, toEstimationScale,
  globals, cdir, cfile, shlib.args,
  userdata, ...,
  .solibs = list(),
  verbose = getOption("verbose",FALSE)
) {

  ep <- paste0("in ",sQuote("pomp"),": ")

  if (missing(t0)) stop(ep,sQuote("t0")," is a required argument",call.=FALSE)

  if (missing(userdata)) userdata <- list()
  added.userdata <- list(...)
  if (length(added.userdata)>0) {
    message(ep,": the following unrecognized argument(s) ",
            "will be stored for use by user-defined functions: ",
            paste(sQuote(names(added.userdata)),collapse=","))
    userdata[names(added.userdata)] <- added.userdata
  }

  if (missing(globals)) globals <- NULL
  if (missing(cdir)) cdir <- NULL
  if (missing(cfile)) cfile <- NULL
  if (missing(shlib.args)) shlib.args <- NULL
  if (missing(PACKAGE)) PACKAGE <- NULL
  PACKAGE <- as.character(PACKAGE)
  globals <- as(globals,"character")

  ## defaults for names of states, parameters, observations, and covariates
  if (missing(statenames)) statenames <- NULL
  if (missing(paramnames)) paramnames <- NULL
  if (missing(obsnames)) obsnames <- NULL
  if (missing(covarnames)) covarnames <- NULL
  if (missing(zeronames)) zeronames <- NULL

  statenames <- as.character(statenames)
  paramnames <- as.character(paramnames)
  zeronames <- as.character(zeronames)
  obsnames <- as.character(obsnames)
  covarnames <- as.character(covarnames)

  ## store the data as double-precision matrix
  storage.mode(data) <- 'double'
  if (length(obsnames)==0) obsnames <- rownames(data)

  ## check for duplicate names
  if (anyDuplicated(statenames)) {
    stop("all ",sQuote("statenames")," must be unique", call.=FALSE)
  }
  if (anyDuplicated(paramnames)) {
    stop("all ",sQuote("paramnames")," must be unique", call.=FALSE)
  }
  if (anyDuplicated(obsnames)) {
    stop("all ",sQuote("paramnames")," must be unique", call.=FALSE)
  }
  if (anyDuplicated(zeronames)) {
    stop("all ",sQuote("zeronames")," must be unique", call.=FALSE)
  }

  ## check the parameters and force them to be double-precision
  if (length(params)>0) {
    if (is.null(names(params)) || !is.numeric(params))
      stop(sQuote("params")," must be a named numeric vector",call.=FALSE)
  }
  storage.mode(params) <- 'double'

  ## check times
  if (!is.numeric(times) || any(is.na(times)) || !all(diff(times)>0))
    stop(sQuote("times")," must be an increasing numeric vector",call.=FALSE)
  storage.mode(times) <- 'double'

  ## check t0
  if (!is.numeric(t0) || length(t0) > 1)
    stop("the zero-time ",sQuote("t0")," must be a single number",call.=FALSE)
  storage.mode(t0) <- 'double'

  ## check and arrange covariates
  if (is.null(covar)) {
    covar <- matrix(data=0,nrow=0,ncol=0)
    tcovar <- numeric(0)
  } else if (is.null(tcovar)) {
    stop("if ",sQuote("covar")," is supplied, ",sQuote("tcovar"),
         " must also be supplied",call.=FALSE)
  } else if (is.data.frame(covar)) {
    if ((is.numeric(tcovar) && (tcovar<1 || tcovar>length(covar))) ||
        (is.character(tcovar) && (!(tcovar%in%names(covar)))) ||
        (!is.numeric(tcovar) && !is.character(tcovar))) {
      stop("if ",sQuote("covar")," is a data frame, ",sQuote("tcovar"),
           " should indicate the time variable",call.=FALSE)
    } else if (is.numeric(tcovar)) {
      tpos <- tcovar
      tcovar <- covar[[tpos]]
      covar <- as.matrix(covar[-tpos])
    } else if (is.character(tcovar)) {
      tpos <- match(tcovar,names(covar))
      tcovar <- covar[[tpos]]
      covar <- as.matrix(covar[-tpos])
    }
  } else {
    covar <- as.matrix(covar)
  }
  if (length(covarnames)==0) covarnames <- as.character(colnames(covar))
  if (!all(covarnames %in% colnames(covar))) {
    missing <- covarnames[!(covarnames%in%colnames(covar))]
    stop("covariate(s) ",paste(sapply(missing,sQuote),collapse=","),
         " are not among the columns of ",sQuote("covar"),call.=FALSE)
  }
  if (anyDuplicated(covarnames)) {
    stop("all ",sQuote("covarnames")," must be unique", call.=FALSE)
  }
  storage.mode(tcovar) <- "double"
  storage.mode(covar) <- "double"

  ## use default initializer?
  default.init <- is.null(initializer) ||
    (is(initializer,"pomp.fun") && initializer@mode == pompfunmode$undef )
  if (default.init) initializer <- pomp.fun(slotname="initializer")

  ## by default, use flat improper prior
  if (is.null(dprior))
    dprior <- pomp.fun(f="_pomp_default_dprior",PACKAGE="pomp")

  ## default rprocess & dprocess
  if (is.null(rprocess))
    rprocess <- function (xstart,times,params,...) stop(sQuote("rprocess")," not specified",call.=FALSE)

  ## handle C snippets
  snips <- list()
  if (is(rprocess,"pompPlugin") && rprocess@csnippet)
    snips <- c(snips,setNames(list(slot(rprocess,rprocess@slotname)@text),rprocess@slotname))
  if (is(dprocess,"Csnippet"))
    snips <- c(snips,dprocess=dprocess@text)
  if (is(skeleton,"Csnippet"))
    snips <- c(snips,skeleton=skeleton@text)
  if (is(rmeasure,"Csnippet"))
    snips <- c(snips,rmeasure=rmeasure@text)
  if (is(dmeasure,"Csnippet"))
    snips <- c(snips,dmeasure=dmeasure@text)
  if (is(rprior,"Csnippet"))
    snips <- c(snips,rprior=rprior@text)
  if (is(dprior,"Csnippet"))
    snips <- c(snips,dprior=dprior@text)
  if (is(fromEstimationScale,"Csnippet"))
    snips <- c(snips,fromEstimationScale=fromEstimationScale@text)
  if (is(toEstimationScale,"Csnippet"))
    snips <- c(snips,toEstimationScale=toEstimationScale@text)
  if (is(initializer,"Csnippet")) {
    if (length(statenames)==0)
      stop(ep,"when ",sQuote("initializer")," is provided as a C snippet, ",
           "you must also provide ",sQuote("statenames"),call.=FALSE)
    snips <- c(snips,initializer=initializer@text)
  }
  if (length(snips)>0) {
    libname <- tryCatch(
      do.call(
        Cbuilder,
        c(
          list(
            name=cfile,
            dir=cdir,
            templates=snippet_templates,
            obsnames=obsnames,
            statenames=statenames,
            paramnames=paramnames,
            covarnames=covarnames,
            globals=globals,
            shlib.args=shlib.args,
            verbose=verbose
          ),
          snips
        )
      ),
      error = function (e) {
        stop("error in building shared-object library from C snippets: ",
             conditionMessage(e),call.=FALSE)
      }
    )
    .solibs <- c(.solibs,list(libname))
    libname <- libname$name
  } else {
    libname <- ''
  }

  ## handle rprocess
  rprocess <- plugin.handler(
    rprocess,
    libname=libname,
    statenames=statenames,
    paramnames=paramnames,
    obsnames=obsnames,
    covarnames=covarnames,
    purpose = sQuote("rprocess")
  )

  ## handle dprocess
  dprocess <- pomp.fun(
    f=dprocess,
    PACKAGE=PACKAGE,
    proto=snippet_templates$dprocess$proto,
    slotname=snippet_templates$dprocess$slotname,
    Cname=snippet_templates$dprocess$Cname,
    libname=libname,
    statenames=statenames,
    paramnames=paramnames,
    obsnames=obsnames,
    covarnames=covarnames
  )

  ## handle initializer
  if (!default.init) {
    initializer <- pomp.fun(
      f=initializer,
      PACKAGE=PACKAGE,
      proto=snippet_templates$initializer$proto,
      slotname=snippet_templates$initializer$slotname,
      Cname=snippet_templates$initializer$Cname,
      libname=libname,
      statenames=statenames,
      paramnames=paramnames,
      obsnames=obsnames,
      covarnames=covarnames
    )
  }

  ## handle skeleton
  if (is.null(skeleton)) skel.type <- "undef"
  skeleton <- pomp.fun(
    f=skeleton,
    PACKAGE=PACKAGE,
    proto=snippet_templates$skeleton$proto,
    slotname=snippet_templates$skeleton$slotname,
    Cname=snippet_templates$skeleton$Cname,
    libname=libname,
    statenames=statenames,
    paramnames=paramnames,
    obsnames=obsnames,
    covarnames=covarnames
  )

  ## type of skeleton (map or vectorfield)
  ## skelmap.delta.t has no meaning in the vectorfield case
  skel.type <- match.arg(skel.type,c("map","vectorfield","undef"))

  ## handle rmeasure
  rmeasure <- pomp.fun(
    f=rmeasure,
    PACKAGE=PACKAGE,
    proto=snippet_templates$rmeasure$proto,
    slotname=snippet_templates$rmeasure$slotname,
    Cname=snippet_templates$rmeasure$Cname,
    libname=libname,
    statenames=statenames,
    paramnames=paramnames,
    obsnames=obsnames,
    covarnames=covarnames
  )

  ## handle dmeasure
  dmeasure <- pomp.fun(
    f=dmeasure,
    PACKAGE=PACKAGE,
    proto=snippet_templates$dmeasure$proto,
    slotname=snippet_templates$dmeasure$slotname,
    Cname=snippet_templates$dmeasure$Cname,
    libname=libname,
    statenames=statenames,
    paramnames=paramnames,
    obsnames=obsnames,
    covarnames=covarnames
  )

  ## handle rprior
  rprior <- pomp.fun(
    f=rprior,
    PACKAGE=PACKAGE,
    proto=snippet_templates$rprior$proto,
    slotname=snippet_templates$rprior$slotname,
    Cname=snippet_templates$rprior$Cname,
    libname=libname,
    statenames=statenames,
    paramnames=paramnames,
    obsnames=obsnames,
    covarnames=covarnames
  )

  ## handle dprior
  dprior <- pomp.fun(
    f=dprior,
    PACKAGE=PACKAGE,
    proto=snippet_templates$dprior$proto,
    slotname=snippet_templates$dprior$slotname,
    Cname=snippet_templates$dprior$Cname,
    libname=libname,
    statenames=statenames,
    paramnames=paramnames,
    obsnames=obsnames,
    covarnames=covarnames
  )

  ## handle parameter transformations
  from.trans <- pomp.fun(
    f=fromEstimationScale,
    PACKAGE=PACKAGE,
    proto=snippet_templates$fromEstimationScale$proto,
    slotname=snippet_templates$fromEstimationScale$slotname,
    Cname=snippet_templates$fromEstimationScale$Cname,
    libname=libname,
    statenames=statenames,
    paramnames=paramnames,
    obsnames=obsnames,
    covarnames=covarnames
  )

  to.trans <- pomp.fun(
    f=toEstimationScale,
    PACKAGE=PACKAGE,
    proto=snippet_templates$toEstimationScale$proto,
    slotname=snippet_templates$toEstimationScale$slotname,
    Cname=snippet_templates$toEstimationScale$Cname,
    libname=libname,
    statenames=statenames,
    paramnames=paramnames,
    obsnames=obsnames,
    covarnames=covarnames
  )

  has.trans <- !is.null(fromEstimationScale) &&
    from.trans@mode != pompfunmode$undef &&
    to.trans@mode != pompfunmode$undef

  ## check to make sure 'covars' is included as an argument where needed
  if (nrow(covar) > 0) {
    if ((skeleton@mode==pompfunmode$Rfun) &&
        !("covars"%in%names(formals(skeleton@R.fun))))
      warning(ep,"a covariate table has been given, yet the ",
              sQuote("skeleton")," function does not have ",
              sQuote("covars")," as a formal argument: see ",
              sQuote("?pomp"),call.=FALSE)
    if ((rmeasure@mode==pompfunmode$Rfun) &&
        !("covars"%in%names(formals(rmeasure@R.fun))))
      warning(ep,"a covariate table has been given, yet the ",
              sQuote("rmeasure")," function does not have ",
              sQuote("covars")," as a formal argument: see ",
              sQuote("?pomp"),call.=FALSE)
    if ((dmeasure@mode==pompfunmode$Rfun) &&
        !("covars"%in%names(formals(dmeasure@R.fun))))
      warning(ep,"a covariate table has been given, yet the ",
              sQuote("dmeasure")," function does not have ",
              sQuote("covars")," as a formal argument: see ",
              sQuote("?pomp"),call.=FALSE)
  }

  if ((length(tcovar)>0)&&((min(tcovar)>t0)||(max(tcovar)<max(times))))
    warning(ep,"the supplied covariate covariate times ",sQuote("tcovar"),
            " do not embrace the data times: covariates may be extrapolated",
            call.=FALSE
    )

  new(
    'pomp',
    rprocess = rprocess,
    dprocess = dprocess,
    dmeasure = dmeasure,
    rmeasure = rmeasure,
    dprior = dprior,
    rprior = rprior,
    skeleton = skeleton,
    skeleton.type = skel.type,
    skelmap.delta.t = skelmap.delta.t,
    data = data,
    times = times,
    t0 = t0,
    default.init = default.init,
    initializer = initializer,
    params = params,
    covar = covar,
    tcovar = tcovar,
    zeronames = zeronames,
    has.trans = has.trans,
    from.trans = from.trans,
    to.trans = to.trans,
    solibs = .solibs,
    userdata = userdata
  )
}

measform2pomp <- function (formulae) {

  ep <- paste0("in ",sQuote("pomp"),": ")

  if (!is.list(formulae))
    formulae <- list(formulae)
  nobs <- length(formulae)
  if (nobs < 1)
    stop(ep,"to use ",sQuote("measurement.model"),
         " you must provide at least one formula",call.=FALSE)
  for (k in seq_len(nobs)) {
    if (!inherits(formulae[[k]],"formula"))
      stop(ep,sQuote("measurement.model")," takes formulae as arguments",call.=FALSE)
  }
  obsnames <- unlist(lapply(formulae,function(x)x[[2L]]))
  distrib <- lapply(formulae,function(x)as.character(x[[3L]][[1L]]))
  ddistrib <- lapply(distrib,function(x)paste0("d",x))
  rdistrib <- lapply(distrib,function(x)paste0("r",x))
  for (k in seq_len(nobs)) {
    tryCatch(
      match.fun(ddistrib[[k]]),
      error = function (e)
        stop(ep,"distribution function ",ddistrib[[k]]," not found",call.=FALSE)
    )
    tryCatch(
      match.fun(rdistrib[[k]]),
      error = function (e)
        stop(ep,"random deviate function ",rdistrib[[k]]," not found",call.=FALSE)
    )
  }
  pred.args <- lapply(formulae,function(x)as.list(x[[3L]][-1L]))
  dcalls <- vector(mode='list',length=nobs)
  rcalls <- vector(mode='list',length=nobs)
  for (k in seq_len(nobs)) {
    dcalls[[k]] <- as.call(
      c(
        list(
          as.name(ddistrib[[k]]),
          x=obsnames[[k]]
        ),
        pred.args[[k]],
        list(
          log=TRUE
        )
      )
    )
    rcalls[[k]] <- as.call(
      c(
        list(
          as.name(rdistrib[[k]]),
          n=1
        ),
        pred.args[[k]]
      )
    )
  }
  list(
    dmeasure = function (y, x, t, params, log, covars, ...) {
      f <- 0
      for (k in seq_len(nobs)) {
        f <- f+eval(
          dcalls[[k]],
          envir=as.list(c(y,x,params,covars,t=t))
        )
      }
      if (log) f else exp(f)
    },
    rmeasure = function (x, t, params, covars, ...) {
      y <- numeric(length=nobs)
      names(y) <- obsnames
      for (k in seq_len(nobs)) {
        y[k] <- eval(
          rcalls[[k]],
          envir=as.list(c(x,params,covars,t=t))
        )
      }
      y
    }
  )
}

vectorfield <- safecall
map <- safecall
