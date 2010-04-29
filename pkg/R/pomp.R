## constructor of the pomp class
default.initializer <- function (params, t0, ...) {
  ivpnames <- grep("\\.0$",names(params),val=TRUE)
  if (length(ivpnames)<1)
    stop("default initializer error: no parameter names ending in ",sQuote(".0")," found: see ",sQuote("pomp")," documentation")
  x <- params[ivpnames]
  names(x) <- sub("\\.0$","",ivpnames)
  x
}

pomp <- function (data, times, t0, ..., rprocess, dprocess,
                  rmeasure, dmeasure, measurement.model,
                  skeleton.map, skeleton.vectorfield, initializer, covar, tcovar,
                  statenames, paramnames, covarnames,
                  PACKAGE) {
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
  
  if (missing(rprocess))
    rprocess <- function(xstart,times,params,...)stop(sQuote("rprocess")," not specified")
  if (missing(dprocess))
    dprocess <- function(x,times,params,log=FALSE,...)stop(sQuote("dprocess")," not specified")

  if (!missing(measurement.model)) {
    if (!(missing(dmeasure)&&missing(rmeasure))) {
      warning("specifying ",sQuote("measurement.model")," overrides specification of ",sQuote("rmeasure")," and ",sQuote("dmeasure"))
    }
    mm <- measform2pomp(measurement.model)
    rmeasure <- mm$rmeasure
    dmeasure <- mm$dmeasure
  }
  
  if (missing(rmeasure))
    rmeasure <- function(x,t,params,covars,...)stop(sQuote("rmeasure")," not specified")
  if (missing(dmeasure))
    dmeasure <- function(y,x,t,params,log=FALSE,covars,...)stop(sQuote("dmeasure")," not specified")
  
  if (missing(skeleton.map)) {
    if (missing(skeleton.vectorfield)) {# skeleton is unspecified
      skeleton.type <- as.character(NA)
      skeleton <- new(
                      "pomp.fun",
                      R.fun=function(x,t,params,covars,...)stop(sQuote("skeleton")," not specified"),
                      use=as.integer(1)
                      )
    } else {                # skeleton is a vectorfield (ordinary differential equation)
      skeleton.type <- "vectorfield"
      if (is.function(skeleton.vectorfield)) {
        if (!all(c('x','t','params','...')%in%names(formals(skeleton.vectorfield))))
          stop(
               sQuote("skeleton.vectorfield"),
               " must be a function of prototype ",
               sQuote("skeleton.vectorfield(x,t,params,...)"),
               call.=TRUE
               )
        skeleton <- new(
                        "pomp.fun",
                        R.fun=skeleton.vectorfield,
                        use=as.integer(1)
                        )
      } else if (is.character(skeleton.vectorfield)) {
        skeleton <- new(
                        "pomp.fun",
                        R.fun=function(x,t,params,...)stop("very bad: skel.fun"),
                        native.fun=skeleton.vectorfield,
                        PACKAGE=PACKAGE,
                        use=as.integer(2)
                        )
      } else {
        stop("pomp error: ",sQuote("skeleton.vectorfield")," must be either a function or the name of a compiled routine")
      }
    }
  } else {
    if (missing(skeleton.vectorfield)) { # skeleton is a map (discrete-time system)
      skeleton.type <- "map"
      if (is.function(skeleton.map)) {
        if (!all(c('x','t','params','...')%in%names(formals(skeleton.map))))
          stop("pomp error: ",sQuote("skeleton.map")," must be a function of prototype ",sQuote("skeleton.map(x,t,params,...)"))
        skeleton <- new(
                        "pomp.fun",
                        R.fun=skeleton.map,
                        use=as.integer(1)
                        )
      } else if (is.character(skeleton.map)) {
        skeleton <- new(
                        "pomp.fun",
                        R.fun=function(x,t,params,...)stop("very bad: skel.fun"),
                        native.fun=skeleton.map,
                        PACKAGE=PACKAGE,
                        use=as.integer(2)
                        )
      } else {
        stop("pomp error: ",sQuote("skeleton.map")," must be either a function or the name of a compiled routine")
      }
    } else { # a dynamical system cannot be both a map and a vectorfield
      stop("pomp error: it is not permitted to specify both ",sQuote("skeleton.map")," and ",sQuote("skeleton.vectorfield"))
    }
  }
  
  if (missing(initializer)) {
    initializer <- default.initializer
  }

  if (missing(PACKAGE)) PACKAGE <- character(0)
  
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

  if (is.function(rmeasure)) {
    if (!all(c('x','t','params','...')%in%names(formals(rmeasure))))
      stop("pomp error: ",sQuote("rmeasure")," must be a function of prototype ",sQuote("rmeasure(x,t,params,...)"))
    rmeasure <- new(
                    "pomp.fun",
                    R.fun=rmeasure,
                    use=as.integer(1)
                    )
  } else if (is.character(rmeasure)) {
    rmeasure <- new(
                    "pomp.fun",
                    R.fun=function(x,t,params,...)stop("very bad: rmeasure.fun"),
                    native.fun=rmeasure,
                    PACKAGE=PACKAGE,
                    use=as.integer(2)
                    )
  } else {
    stop("pomp error: ",sQuote("rmeasure")," must be either a function or the name of a compiled routine")
  }
  
  if (is.function(dmeasure)) {
    if (!all(c('y','x','t','params','log','...')%in%names(formals(dmeasure))))
      stop("pomp error: ",sQuote("dmeasure")," must be a function of prototype ",sQuote("dmeasure(y,x,t,params,log,...)"))
    dmeasure <- new(
                    "pomp.fun",
                    R.fun=dmeasure,
                    use=as.integer(1)
                    )
  } else if (is.character(dmeasure)) {
    dmeasure <- new(
                    "pomp.fun",
                    R.fun=function(y,x,t,params,log,...)stop("very bad: rmeasure.fun"),
                    native.fun=dmeasure,
                    PACKAGE=PACKAGE,
                    use=as.integer(2)
                    )
  } else {
    stop("pomp error: ",sQuote("dmeasure")," must be either a function or the name of a compiled routine")
  }
  
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
  
  if (missing(statenames)) statenames <- character(0)
  if (missing(paramnames)) paramnames <- character(0)
  if (missing(covarnames)) covarnames <- character(0)
  
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
      warning("a covariate table has been given, yet the ",sQuote("skeleton")," function does not have ",sQuote("covars")," as a formal argument")
    if (
        (rmeasure@use==1)
        &&!("covars"%in%names(formals(rmeasure@R.fun)))
        )
      warning("a covariate table has been given, yet the ",sQuote("rmeasure")," function does not have ",sQuote("covars")," as a formal argument")
    if (
        (dmeasure@use==1)
        &&!("covars"%in%names(formals(dmeasure@R.fun)))
        )
      warning("a covariate table has been given, yet the ",sQuote("dmeasure")," function does not have ",sQuote("covars")," as a formal argument")
  }

  new(
      'pomp',
      rprocess = rprocess,
      dprocess = dprocess,
      dmeasure = dmeasure,
      rmeasure = rmeasure,
      skeleton.type = skeleton.type,
      skeleton = skeleton,
      data = data,
      times = times,
      t0 = t0,
      initializer = initializer,
      covar = covar,
      tcovar = tcovar,
      statenames = statenames,
      paramnames = paramnames,
      covarnames = covarnames,
      PACKAGE = PACKAGE,
      userdata = list(...)
      )
}

measform2pomp <- function (formulae) {
  if (!is.list(formulae))
    formulae <- list(formulae)
  nobs <- length(formulae)
  if (nobs < 1)
    stop("pomp error: to use ",sQuote("measurement.model")," you must provide at least one formula")
  for (k in seq_len(nobs)) {
    if (!inherits(formulae[[k]],"formula"))
      stop("pomp error: ",sQuote("measurement.model")," takes formulae as arguments")
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

