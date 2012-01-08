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
                        skeleton.type = 'character',
                        skeleton = 'pomp.fun',
                        skelmap.delta.t = 'numeric',
                        initializer = 'function',
                        states = 'array',
                        params = 'numeric',
                        covar = 'matrix',
                        tcovar = 'numeric',
                        obsnames = 'character',
                        statenames = 'character',
                        paramnames = 'character',
                        covarnames = 'character',
                        zeronames = 'character',
                        par.trans = 'function',
                        par.untrans = 'function',
                        PACKAGE = 'character',
                        userdata = 'list'
                        )
         )

## this is the initial-condition setting function that is used by default
## it simply finds all parameters in the vector 'params' that have a name ending in '.0'
## and returns a vector with their values with names stripped of '.0'
default.initializer <- function (params, t0, ...) {
  ivpnames <- grep("\\.0$",names(params),value=TRUE)
  if (length(ivpnames)<1)
    stop("default initializer error: no parameter names ending in ",sQuote(".0")," found: see ",sQuote("pomp")," documentation")
  x <- params[ivpnames]
  names(x) <- sub("\\.0$","",ivpnames)
  x
}

## as of version 0.37-1 'pomp' is a generic function
setGeneric("pomp",function(data,...)standardGeneric("pomp"))

## basic constructor of the pomp class
pomp.constructor <- function (data, times, t0, ..., rprocess, dprocess,
                              rmeasure, dmeasure, measurement.model,
                              skeleton = NULL, skeleton.type = c("map","vectorfield"),
                              skelmap.delta.t = 1,
                              initializer, covar, tcovar,
                              obsnames, statenames, paramnames, covarnames, zeronames,
                              PACKAGE, parameter.transform, parameter.inv.transform) {

  ## check the data
  if (is.data.frame(data)) {
    if (!is.character(times) || length(times)!=1 || !(times%in%names(data)))
      stop("pomp error: ",sQuote("times")," must be the name of a column of ",sQuote("data"),call.=TRUE)
    tmnm <- times
    times <- data[[tmnm]]
    data <- do.call(rbind,lapply(data[!(names(data)%in%tmnm)],as.numeric))
  }
  if (!is.numeric(data))
    stop("pomp error: ",sQuote("data")," must be numeric",call.=TRUE)
  if (!is.array(data))
    data <- array(data,dim=c(1,length(data)),dimnames=list("data",NULL))
  storage.mode(data) <- 'double'
  
  ## check times
  if (!is.numeric(times) || !all(diff(times)>0))
    stop("pomp error: ",sQuote("times")," must be an increasing numeric vector",call.=TRUE)
  if (length(times)!=ncol(data))
    stop("pomp error: the length of ",sQuote("times")," does not equal the number of columns in ",sQuote("data"),call.=TRUE)
  storage.mode(times) <- 'double'
  
  ## check t0
  if (!is.numeric(t0) || length(t0) > 1)
    stop("pomp error: the zero-time ",sQuote("t0")," must be a single number",call.=TRUE)
  if (t0 > times[1])
    stop("pomp error: the zero-time ",sQuote("t0")," must occur no later than the first observation",call.=TRUE)
  storage.mode(t0) <- 'double'
  
  if (missing(PACKAGE)) PACKAGE <- ""
  
  if (missing(rprocess))
    rprocess <- function(xstart,times,params,...)stop(sQuote("rprocess")," not specified")
  if (missing(dprocess))
    dprocess <- function(x,times,params,log=FALSE,...)stop(sQuote("dprocess")," not specified")

  if (!missing(measurement.model)) {
    if (!(missing(dmeasure)&&missing(rmeasure))) {
      warning(
              "specifying ",sQuote("measurement.model"),
              " overrides specification of ",sQuote("rmeasure")," and ",sQuote("dmeasure")
              )
    }
    mm <- measform2pomp(measurement.model)
    rmeasure <- mm$rmeasure
    dmeasure <- mm$dmeasure
  }
  
  if (missing(rmeasure))
    rmeasure <- function(x,t,params,covars,...)stop(sQuote("rmeasure")," not specified")
  if (missing(dmeasure))
    dmeasure <- function(y,x,t,params,log,covars,...)stop(sQuote("dmeasure")," not specified")
  
  skeleton.type <- match.arg(skeleton.type)
  skelmap.delta.t <- as.numeric(skelmap.delta.t)
  if (skelmap.delta.t <= 0)
    stop(sQuote("skelmap.delta.t")," must be positive")

  if (is.null(skeleton)) {
    skeleton <- pomp.fun(f=function(x,t,params,covars,...)stop(sQuote("skeleton")," not specified"))
  } else {
    skeleton <- pomp.fun(f=skeleton,PACKAGE=PACKAGE,proto=quote(skeleton(x,t,params,...)))
  }
  
  if (missing(initializer)) {
    initializer <- default.initializer
  }

  if (!is.function(rprocess))
    stop(
         "pomp error: ",sQuote("rprocess")," must be a function",
         call.=TRUE
         )
  if (!is.function(dprocess))
    stop(
         "pomp error: ",sQuote("dprocess")," must be a function",
         call.=TRUE
         )

  rmeasure <- pomp.fun(f=rmeasure,PACKAGE=PACKAGE,proto=quote(rmeasure(x,t,params,...)))
  dmeasure <- pomp.fun(f=dmeasure,PACKAGE=PACKAGE,proto=quote(dmeasure(y,x,t,params,log,...)))
  
  if (!is.function(initializer))
    stop(
         "pomp error: ",sQuote("initializer")," must be a function",
         call.=TRUE
         )
  
  if (!all(c('xstart','times','params','...')%in%names(formals(rprocess))))
    stop(
         "pomp error: ",sQuote("rprocess")," must be a function of prototype ",sQuote("rprocess(xstart,times,params,...)"),
         call.=TRUE
         )
  if (!all(c('x','times','params','log','...')%in%names(formals(dprocess))))
    stop(
         "pomp error: ",sQuote("dprocess")," must be a function of prototype ",sQuote("dprocess(x,times,params,log,...)"),
         call.=TRUE
         )
  if (!all(c('params','t0','...')%in%names(formals(initializer))))
    stop(
         "pomp error: ",sQuote("initializer")," must be a function of prototype ",sQuote("initializer(params,t0,...)"),
         call.=TRUE
         )
  
  if (missing(obsnames)) obsnames <- character(0)
  if (missing(statenames)) statenames <- character(0)
  if (missing(paramnames)) paramnames <- character(0)
  if (missing(covarnames)) covarnames <- character(0)
  if (missing(zeronames)) zeronames <- character(0)
  
  if (missing(covar)) {
    covar <- matrix(data=0,nrow=0,ncol=0)
    tcovar <- numeric(0)
    covarnames <- character(0)
  } else if (missing(tcovar)) {
    stop("pomp error: if ",sQuote("covar")," is supplied, ",sQuote("tcovar")," must be supplied also")
  } else if (is.data.frame(covar)) {
    if (
        (
         is.numeric(tcovar)&&
         (
          (tcovar<1)||
          (tcovar>length(covar))
          )
         )||
        (
         is.character(tcovar)&&
         !(tcovar%in%names(covar))
         )
        ) {
      stop("pomp error: if ",sQuote("covar")," is a data frame, ",sQuote("tcovar")," should indicate the time variable")
    } else {
      tpos <- match(tcovar,names(covar))
      tcovar <- covar[[tpos]]
      covar <- as.matrix(covar[-tpos])
    }
  } else {
    covar <- as.matrix(covar)
  }
  
  storage.mode(tcovar) <- "double"
  storage.mode(covar) <- "double"

  if (length(tcovar)!=nrow(covar)) {
    stop("pomp error: the length of ",sQuote("tcovar")," should match the number of rows of ",sQuote("covar"))
  } else if (!all(covarnames%in%colnames(covar))) {
    missing <- covarnames[!(covarnames%in%colnames(covar))]
    stop("pomp error: covariate(s) ",paste(missing,collapse=",")," are not found among the columns of ",sQuote("covar"))
  } else if (!is.numeric(tcovar)) {
    stop("pomp error: ",sQuote("tcovar")," must either be a numeric vector or must name a numeric vector in the data frame ",sQuote("covar"))
  }

  if (nrow(covar)>0) {
    if (
        (skeleton@use==1)
        &&!("covars"%in%names(formals(skeleton@R.fun)))
        )
      warning("a covariate table has been given, yet the ",sQuote("skeleton")," function does not have ",sQuote("covars")," as a formal argument",call.=FALSE)
    if (
        (rmeasure@use==1)
        &&!("covars"%in%names(formals(rmeasure@R.fun)))
        )
      warning("a covariate table has been given, yet the ",sQuote("rmeasure")," function does not have ",sQuote("covars")," as a formal argument",call.=FALSE)
    if (
        (dmeasure@use==1)
        &&!("covars"%in%names(formals(dmeasure@R.fun)))
        )
      warning("a covariate table has been given, yet the ",sQuote("dmeasure")," function does not have ",sQuote("covars")," as a formal argument",call.=FALSE)
  }

  if ((length(tcovar)>0)&&((min(tcovar)>t0)||(max(tcovar)<max(times)))) 
    warning(
            "the supplied covariate covariate times ",sQuote("tcovar"),
            " do not embrace the data times: covariates may be extrapolated"
            )

  if (missing(parameter.transform)) {
    if (missing(parameter.inv.transform)) {
      has.trans <- FALSE
    } else {
      stop("pomp error: if ",sQuote("parameter.inv.transform")," is supplied, then " ,
           sQuote("parameter.transform")," must also be supplied")
    }
  } else {
    if (missing(parameter.inv.transform)) {
      stop("pomp error: if ",sQuote("parameter.transform")," is supplied, then " ,
           sQuote("parameter.inv.transform")," must also be supplied")
    } else {
      has.trans <- TRUE
    }
  }
  if (has.trans) {
    par.trans <- match.fun(parameter.transform)
    par.untrans <- match.fun(parameter.inv.transform)
    if (!all(c('params','...')%in%names(formals(par.trans))))
      stop(
           "pomp error: ",sQuote("parameter.transform")," must be a function of prototype ",
           sQuote("parameter.transform(params,...)"),
           call.=TRUE
           )
    if (!all(c('params','...')%in%names(formals(par.untrans))))
      stop(
           "pomp error: ",sQuote("parameter.inv.transform")," must be a function of prototype ",
           sQuote("parameter.inv.transform(params,...)"),
           call.=TRUE
           )
  } else {
    par.untrans <- par.trans <- function(params, ...) params
  }
  

  new(
      'pomp',
      rprocess = rprocess,
      dprocess = dprocess,
      dmeasure = dmeasure,
      rmeasure = rmeasure,
      skeleton = skeleton,
      skeleton.type = skeleton.type,
      skelmap.delta.t=skelmap.delta.t,
      data = data,
      times = times,
      t0 = t0,
      initializer = initializer,
      covar = covar,
      tcovar = tcovar,
      obsnames = obsnames,
      statenames = statenames,
      paramnames = paramnames,
      covarnames = covarnames,
      zeronames = zeronames,
      par.trans = par.trans,
      par.untrans = par.untrans,
      PACKAGE = PACKAGE,
      userdata = list(...)
      )
}

measform2pomp <- function (formulae) {
  if (!is.list(formulae))
    formulae <- list(formulae)
  nobs <- length(formulae)
  if (nobs < 1)
    stop("pomp error: to use ",sQuote("measurement.model")," you must provide at least one formula",call.=FALSE)
  for (k in seq_len(nobs)) {
    if (!inherits(formulae[[k]],"formula"))
      stop("pomp error: ",sQuote("measurement.model")," takes formulae as arguments",call.=FALSE)
  }
  obsnames <- unlist(lapply(formulae,function(x)x[[2]]))
  distrib <- lapply(formulae,function(x)as.character(x[[3]][[1]]))
  ddistrib <- lapply(distrib,function(x)paste("d",x,sep=''))
  rdistrib <- lapply(distrib,function(x)paste("r",x,sep=''))
  for (k in seq_len(nobs)) {
    res <- try(
               match.fun(ddistrib[[k]]),
               silent=TRUE
               )
    if (inherits(res,'try-error'))
      stop("pomp error: distribution function ",ddistrib[[k]]," not found")
    res <- try(
               match.fun(rdistrib[[k]]),
               silent=TRUE
               )
    if (inherits(res,'try-error'))
      stop("pomp error: random deviate function ",rdistrib[[k]]," not found")
  }
  pred.args <- lapply(formulae,function(x)as.list(x[[3]][-1]))
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

setMethod(
          "pomp",
          signature(data="data.frame"),
          function (data, times, t0, ..., rprocess, dprocess,
                    rmeasure, dmeasure, measurement.model,
                    skeleton = NULL, skeleton.type = c("map","vectorfield"),
                    skelmap.delta.t = 1,
                    initializer, covar, tcovar,
                    obsnames, statenames, paramnames, covarnames, zeronames,
                    PACKAGE, parameter.transform, parameter.inv.transform) {
            pomp.constructor(
                             data=data,
                             times=times,
                             t0=t0,
                             rprocess=rprocess,
                             dprocess=dprocess,
                             rmeasure=rmeasure,
                             dmeasure=dmeasure,
                             measurement.model=measurement.model,
                             skeleton=skeleton,
                             skeleton.type=skeleton.type,
                             skelmap.delta.t=skelmap.delta.t,
                             initializer=initializer,
                             covar=covar,
                             tcovar=tcovar,
                             obsnames=obsnames,
                             statenames=statenames,
                             paramnames=paramnames,
                             covarnames=covarnames,
                             zeronames=zeronames,
                             PACKAGE=PACKAGE,
                             parameter.transform=parameter.transform,
                             parameter.inv.transform=parameter.inv.transform,
                             ...
                             )
          }
          )

setMethod(
          "pomp",
          signature(data="matrix"),
          function (data, times, t0, ..., rprocess, dprocess,
                    rmeasure, dmeasure, measurement.model,
                    skeleton = NULL, skeleton.type = c("map","vectorfield"),
                    skelmap.delta.t = 1,
                    initializer, covar, tcovar,
                    obsnames, statenames, paramnames, covarnames, zeronames,
                    PACKAGE, parameter.transform, parameter.inv.transform) {
            pomp.constructor(
                             data=data,
                             times=times,
                             t0=t0,
                             rprocess=rprocess,
                             dprocess=dprocess,
                             rmeasure=rmeasure,
                             dmeasure=dmeasure,
                             measurement.model=measurement.model,
                             skeleton=skeleton,
                             skeleton.type=skeleton.type,
                             skelmap.delta.t=skelmap.delta.t,
                             initializer=initializer,
                             covar=covar,
                             tcovar=tcovar,
                             obsnames=obsnames,
                             statenames=statenames,
                             paramnames=paramnames,
                             covarnames=covarnames,
                             zeronames=zeronames,
                             PACKAGE=PACKAGE,
                             parameter.transform=parameter.transform,
                             parameter.inv.transform=parameter.inv.transform,
                             ...
                             )
          }
          )


setMethod(
          "pomp",
          signature(data="numeric"),
          function (data, times, t0, ..., rprocess, dprocess,
                    rmeasure, dmeasure, measurement.model,
                    skeleton = NULL, skeleton.type = c("map","vectorfield"),
                    skelmap.delta.t = 1,
                    initializer, covar, tcovar,
                    obsnames, statenames, paramnames, covarnames, zeronames,
                    PACKAGE, parameter.transform, parameter.inv.transform) {
            pomp.constructor(
                             data=matrix(data,nrow=1,ncol=length(data)),
                             times=times,
                             t0=t0,
                             rprocess=rprocess,
                             dprocess=dprocess,
                             rmeasure=rmeasure,
                             dmeasure=dmeasure,
                             measurement.model=measurement.model,
                             skeleton=skeleton,
                             skeleton.type=skeleton.type,
                             skelmap.delta.t=skelmap.delta.t,
                             initializer=initializer,
                             covar=covar,
                             tcovar=tcovar,
                             obsnames=obsnames,
                             statenames=statenames,
                             paramnames=paramnames,
                             covarnames=covarnames,
                             zeronames=zeronames,
                             PACKAGE=PACKAGE,
                             parameter.transform=parameter.transform,
                             parameter.inv.transform=parameter.inv.transform,
                             ...
                             )
          }
          )

setMethod(
          "pomp",
          signature(data="pomp"),
          function (data, times, t0, ..., rprocess, dprocess,
                    rmeasure, dmeasure, measurement.model,
                    skeleton, skeleton.type, skelmap.delta.t,
                    initializer, covar, tcovar,
                    obsnames, statenames, paramnames, covarnames, zeronames,
                    PACKAGE, parameter.transform, parameter.inv.transform) {
            mmg <- !missing(measurement.model)
            dmg <- !missing(dmeasure)
            rmg <- !missing(rmeasure)
            if (missing(times)) times <- data@times
            if (missing(t0)) t0 <- data@t0
            if (mmg) {
              if (dmg||rmg)
                warning(
                        "specifying ",sQuote("measurement.model"),
                        " overrides specification of ",sQuote("rmeasure")," and ",sQuote("dmeasure")
                        )
              mm <- measform2pomp(measurement.model)
              rmeasure <- mm$rmeasure
              dmeasure <- mm$dmeasure
            } else {
              if (!rmg) rmeasure <- data@rmeasure
              if (!dmg) dmeasure <- data@dmeasure
            }
            if (missing(rprocess)) rprocess <- data@rprocess
            if (missing(dprocess)) dprocess <- data@dprocess
            if (missing(initializer)) initializer <- data@initializer
            if (missing(covar)) covar <- data@covar
            if (missing(tcovar)) tcovar <- data@tcovar
            if (missing(obsnames)) obsnames <- data@obsnames
            if (missing(statenames)) statenames <- data@statenames
            if (missing(paramnames)) paramnames <- data@paramnames
            if (missing(covarnames)) covarnames <- data@covarnames
            if (missing(zeronames)) zeronames <- data@zeronames
            if (missing(PACKAGE)) PACKAGE <- data@PACKAGE
            if (missing(skeleton.type)) skeleton.type <- data@skeleton.type
            if (missing(skeleton)) skeleton <- data@skeleton
            if (missing(skelmap.delta.t)) skelmap.delta.t <- data@skelmap.delta.t

            if (missing(parameter.transform)) {
              if (missing(parameter.inv.transform)) {
                par.trans <- data@par.trans
                par.untrans <- data@par.untrans
              } else {
                stop("pomp error: if ",sQuote("parameter.inv.transform")," is supplied, then " ,
                     sQuote("parameter.transform")," must also be supplied")
              }
            } else {
              if (missing(parameter.inv.transform)) {
                stop("pomp error: if ",sQuote("parameter.transform")," is supplied, then " ,
                     sQuote("parameter.inv.transform")," must also be supplied")
              } else {
                par.trans <- match.fun(parameter.transform)
                par.untrans <- match.fun(parameter.inv.transform)
              }
            }
            
            userdata <- data@userdata
            added.userdata <- list(...)
            userdata[names(added.userdata)] <- added.userdata
            do.call(
                    pomp.constructor,
                    c(
                      list(
                           data=data@data,
                           times=times,
                           t0=t0,
                           rprocess=rprocess,
                           dprocess=dprocess,
                           rmeasure=rmeasure,
                           dmeasure=dmeasure,
                           skeleton=skeleton,
                           skeleton.type=skeleton.type,
                           skelmap.delta.t=skelmap.delta.t,
                           initializer=initializer,
                           covar=covar,
                           tcovar=tcovar,
                           obsnames=obsnames,
                           statenames=statenames,
                           paramnames=paramnames,
                           covarnames=covarnames,
                           zeronames=zeronames,
                           PACKAGE=PACKAGE,
                           parameter.transform=par.trans,
                           parameter.inv.transform=par.untrans
                           ),
                      userdata
                      )
                    )
          }
          )
