## This file defines 'pomp', the basic constructor of the pomp class

pomp <- function (data, times, t0, ..., rprocess, dprocess,
                  rmeasure, dmeasure, measurement.model,
                  skeleton, initializer, rprior, dprior, params, covar, tcovar,
                  obsnames, statenames, paramnames, covarnames, zeronames,
                  PACKAGE, fromEstimationScale, toEstimationScale,
                  globals, cdir, cfile, shlib.args) {

  ep <- paste0("in ",sQuote("pomp"),": ")  # error prefix

  if (missing(data))
    stop(ep,sQuote("data")," is a required argument",call.=FALSE)

  ## return as quickly as possible if no work is to be done
  if (nargs()==1 && is(data,"pomp")) return(data)

  ## if one transformation is supplied, then both must be
  fes <- !missing(fromEstimationScale)
  tes <- !missing(toEstimationScale)
  if (xor(fes,tes))
    stop(ep,"if one of ",sQuote("toEstimationScale"),", ",
         sQuote("fromEstimationScale")," is supplied, then so must the other",
         call.=FALSE)

  ## if 'measurement model' is specified as a formula, this overrides
  ## specification of 'rmeasure' or 'dmeasure'
  mmg <- !missing(measurement.model)
  dmg <- !(missing(dmeasure) || is.null(dmeasure))
  rmg <- !(missing(rmeasure) || is.null(rmeasure))
  if (mmg) {
    if (dmg||rmg)
      warning(ep,"specifying ",sQuote("measurement.model"),
              " overrides specification of ",
              sQuote("rmeasure")," and ",sQuote("dmeasure"),
              call.=FALSE)
    mm <- measform2pomp(measurement.model)
    rmeasure <- mm$rmeasure
    dmeasure <- mm$dmeasure
    dmg <- TRUE
    rmg <- TRUE
  }

  ## the deterministic skeleton involves 'skeleton', 'skeleton.type', and
  ## 'skelmap.delta.t'
  skel.type <- "undef"
  skel.delta.t <- 1

  has.skel <- !missing(skeleton)
  if (has.skel) {
    skeleton <- substitute(skeleton)
    flist <- list(
      map=function (f, delta.t = 1) {
        skel.type <<- "map"
        if (delta.t <= 0)
          stop("in ",sQuote("map"),", ",sQuote("delta.t"),
               " must be positive",call.=FALSE)
        skel.delta.t <<- as.numeric(delta.t)
        f
      },
      vectorfield=function (f) {
        skel.type <<- "vectorfield"
        f
      }
    )
    skeleton <- eval(skeleton,envir=flist,enclos=parent.frame())
  }

  ## by default, use flat improper prior
  if (missing(dprior))
    dprior <- pomp.fun(f="_pomp_default_dprior",PACKAGE="pomp")

  if (missing(globals)) globals <- NULL
  if (missing(cdir)) cdir <- NULL
  if (missing(cfile)) cfile <- NULL
  if (missing(shlib.args)) shlib.args <- NULL
  if (missing(PACKAGE)) PACKAGE <- NULL

  ## defaults for names of states, parameters, observations, and covariates
  if (missing(statenames)) statenames <- NULL
  if (missing(paramnames)) paramnames <- NULL
  if (missing(obsnames)) obsnames <- NULL
  if (missing(covarnames)) covarnames <- NULL

  if (is(data,"pomp")) {
    ## 'data' is a pomp object:
    ## extract missing arguments from it

    if (missing(times)) times <- data@times
    if (missing(t0)) t0 <- data@t0

    if (missing(rmeasure)) rmeasure <- data@rmeasure
    if (missing(dmeasure)) dmeasure <- data@dmeasure

    if (missing(rprocess)) rprocess <- data@rprocess
    if (missing(dprocess)) dprocess <- data@dprocess
    if (missing(rprior)) rprior <- data@rprior
    if (missing(initializer)) initializer <- data@initializer

    if (missing(params)) params <- data@params
    if (missing(covar)) covar <- data@covar
    if (missing(tcovar)) tcovar <- data@tcovar
    if (missing(zeronames)) zeronames <- data@zeronames

    ## the deterministic skeleton involves 'skeleton', 'skeleton.type', and
    ## 'skelmap.delta.t'
    if (!has.skel) {
      skel.type <- data@skeleton.type
      skel.delta.t <- data@skelmap.delta.t
      skeleton <- data@skeleton
    }

    if (fes && tes) {
      from.trans <- fromEstimationScale
      to.trans <- toEstimationScale
    } else {
      from.trans <- data@from.trans
      to.trans <- data@to.trans
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
        skel.delta.t=skel.delta.t,
        initializer=initializer,
        covar=covar,
        tcovar=tcovar,
        obsnames=obsnames,
        statenames=statenames,
        paramnames=paramnames,
        covarnames=covarnames,
        zeronames=zeronames,
        PACKAGE=PACKAGE,
        fromEstimationScale=from.trans,
        toEstimationScale=to.trans,
        params=params,
        globals=globals,
        cdir=cdir,
        cfile=cfile,
        shlib.args=shlib.args,
        .solibs=data@solibs,
        userdata=data@userdata,
        ...
      ),
      error = function (e) {
        stop(ep,conditionMessage(e),call.=FALSE)
      }
    )
  } else if (is.data.frame(data)) {
    ## 'data' is a data frame:
    ## construct a pomp object de novo

    ## preliminary error checking
    if (missing(times)) stop(sQuote("times")," is a required argument",call.=FALSE)
    if (missing(t0)) stop(sQuote("t0")," is a required argument",call.=FALSE)
    if ((is.numeric(times) && (times<1 || times>ncol(data) ||times!=as.integer(times))) ||
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

    if (missing(initializer)) initializer <- NULL
    if (missing(rprocess)) rprocess <- NULL
    if (missing(dprocess)) dprocess <- NULL
    if (missing(rprior)) rprior <- NULL

    if (missing(params)) params <- numeric(0)
    if (missing(covar)) covar <- NULL
    if (missing(tcovar)) tcovar <- NULL
    if (missing(zeronames)) zeronames <- NULL

    if (!has.skel) skeleton <- NULL
    if (!rmg) rmeasure <- NULL
    if (!dmg) dmeasure <- NULL
    if (fes && tes) {
      from.trans <- fromEstimationScale
      to.trans <- toEstimationScale
    } else {
      from.trans <- NULL
      to.trans <- NULL
    }

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
        skel.delta.t=skel.delta.t,
        initializer=initializer,
        params=params,
        covar=covar,
        tcovar=tcovar,
        obsnames=obsnames,
        statenames=statenames,
        paramnames=paramnames,
        covarnames=covarnames,
        zeronames=zeronames,
        PACKAGE=PACKAGE,
        fromEstimationScale=from.trans,
        toEstimationScale=to.trans,
        globals=globals,
        cdir=cdir,
        cfile=cfile,
        shlib.args=shlib.args,
        ...
      ),
      error = function (e) {
        stop(ep,conditionMessage(e),call.=FALSE)
      }
    )
  } else {
    stop(ep,sQuote("data"),
         " must be a data frame or an object of class ",sQuote("pomp"),
         call.=FALSE)
  }
}

pomp.internal <- function (data, times, t0, rprocess, dprocess,
                           rmeasure, dmeasure,
                           skeleton, skel.type, skel.delta.t,
                           initializer, rprior, dprior,
                           params, covar, tcovar,
                           obsnames, statenames, paramnames, covarnames,
                           zeronames, PACKAGE,
                           fromEstimationScale, toEstimationScale,
                           globals, cdir, cfile, shlib.args,
                           userdata, ...,
                           .solibs = list(),
                           verbose = getOption("verbose",FALSE)) {

  ep <- paste0("in ",sQuote("pomp"),": ")

  if (missing(userdata)) userdata <- list()
  added.userdata <- list(...)
  if (length(added.userdata)>0) {
    message(ep,": the following unrecognized argument(s) ",
            "will be stored for use by user-defined functions: ",
            paste(sQuote(names(added.userdata)),collapse=","))
    userdata[names(added.userdata)] <- added.userdata
  }

  ## name of shared object library
  PACKAGE <- as.character(PACKAGE)

  if (!is(globals,"Csnippet")) globals <- as.character(globals)

  statenames <- as.character(statenames)
  paramnames <- as.character(paramnames)
  zeronames <- as.character(zeronames)
  obsnames <- as.character(obsnames)
  covarnames <- as.character(covarnames)

  ## check for duplicate names
  if (anyDuplicated(statenames)) {
    stop("all ",sQuote("statenames")," must be unique", call.=FALSE)
  }
  if (anyDuplicated(paramnames)) {
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

  ## store the data as double-precision matrix
  storage.mode(data) <- 'double'
  if (length(obsnames)==0) obsnames <- rownames(data)

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
  storage.mode(tcovar) <- "double"
  storage.mode(covar) <- "double"

  ## use default initializer?
  default.init <- is.null(initializer) ||
    ( is(initializer,"pomp.fun") && initializer@mode == pompfunmode$undef )
  if (default.init) initializer <- pomp.fun(slotname="initializer")

  ## default rprocess & dprocess
  if (is.null(rprocess))
    rprocess <- function (xstart,times,params,...) stop(sQuote("rprocess")," not specified",call.=FALSE)
  if (is.null(dprocess))
    dprocess <- function (x,times,params,log=FALSE,...) stop(sQuote("dprocess")," not specified",call.=FALSE)

  ## handle C snippets
  snips <- list()
  if (is(rprocess,"pompPlugin") && rprocess@csnippet)
    snips <- c(snips,setNames(list(slot(rprocess,rprocess@slotname)@text),rprocess@slotname))
  if (is(dprocess,"pompPlugin") && dprocess@csnippet)
    snips <- c(snips,setNames(list(slot(dprocess,dprocess@slotname)@text),dprocess@slotname))
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
  dprocess <- plugin.handler(
    dprocess,
    libname=libname,
    statenames=statenames,
    paramnames=paramnames,
    obsnames=obsnames,
    covarnames=covarnames,
    purpose = sQuote("dprocess")
  )

  ## handle initializer
  if (!default.init) {
    initializer <- pomp.fun(
      f=initializer,
      PACKAGE=PACKAGE,
      proto=snippet_templates$initializer$proto,
      slotname="initializer",
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
    slotname="skeleton",
    Cname=snippet_templates$skeleton$Cname,
    libname=libname,
    statenames=statenames,
    paramnames=paramnames,
    obsnames=obsnames,
    covarnames=covarnames
  )

  ## type of skeleton (map or vectorfield)
  ## skel.delta.t has no meaning in the vectorfield case
  skel.type <- match.arg(skel.type,c("map","vectorfield","undef"))
  skel.delta.t <- as.numeric(skel.delta.t)
  if (skel.delta.t <= 0)
    stop("skeleton ",sQuote("delta.t")," must be positive",call.=FALSE) #nocov

  ## handle rmeasure
  rmeasure <- pomp.fun(
    f=rmeasure,
    PACKAGE=PACKAGE,
    proto=snippet_templates$rmeasure$proto,
    slotname="rmeasure",
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
    slotname="dmeasure",
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
    proto=quote(rprior(params,...)),
    slotname="rprior",
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
    slotname="dprior",
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
    slotname="fromEstimationScale",
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
    slotname="toEstimationScale",
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
    skelmap.delta.t = skel.delta.t,
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
