## constructor of the pomp class
pomp <- function (data, times, t0, ..., rprocess, dprocess,
                  rmeasure, dmeasure, measurement.model,
                  skeleton, initializer, covar, tcovar,
                  statenames, paramnames, covarnames,
                  PACKAGE) {
  ## check the data
  if (is.data.frame(data)) {
    if (!is.character(times) || length(times)!=1 || !(times%in%names(data)))
      stop("pomp error: 'times' must be the name of a column of 'data'")
    tmnm <- times
    times <- data[[tmnm]]
    data <- do.call(rbind,lapply(data[!(names(data)%in%tmnm)],as.numeric))
  }
  if (!is.numeric(data))
    stop("pomp error: 'data' must be numeric")
  if (!is.array(data))
    data <- array(data,dim=c(1,length(data)),dimnames=list("data",NULL))
  storage.mode(data) <- 'double'
  
  ## check times
  if (!is.numeric(times) || !all(diff(times)>0))
    stop("pomp error: 'times' must be an increasing numeric vector")
  if (length(times)!=ncol(data))
    stop("pomp error: the length of 'times' does not equal the number of columns in 'data'")
  storage.mode(times) <- 'double'
  
  ## check t0
  if (!is.numeric(t0) || length(t0) > 1)
    stop("pomp error: the zero-time 't0' must be a single number")
  if (t0 > times[1])
    stop("pomp error: the zero-time 't0' must occur no later than the first observation")
  storage.mode(t0) <- 'double'
  
  if (missing(rprocess))
    rprocess <- function(xstart,times,params,...)stop("'rprocess' not specified")
  if (missing(dprocess))
    dprocess <- function(x,times,params,log=FALSE,...)stop("'dprocess' not specified")

  if (!missing(measurement.model)) {
    if (!(missing(dmeasure)&&missing(rmeasure))) {
      warning("specifying 'measurement.model' overrides specification of 'rmeasure' and 'dmeasure'")
    }
    mm <- measform2pomp(measurement.model)
    rmeasure <- mm$rmeasure
    dmeasure <- mm$dmeasure
  }
  
  if (missing(rmeasure))
    rmeasure <- function(x,t,params,covars,...)stop("'rmeasure' not specified")
  if (missing(dmeasure))
    dmeasure <- function(y,x,t,params,log=FALSE,covars,...)stop("'dmeasure' not specified")
  if (missing(skeleton))
    skeleton <- function(x,t,params,covars,...)stop("'skeleton' not specified")

  if (missing(initializer)) {
    initializer <- function (params, t0, ...) {
      ivpnames <- grep("\\.0$",names(params),val=TRUE)
      if (length(ivpnames)<1)
        stop("default initializer error: no parameter names ending in '.0' found: see 'pomp' documentation")
      x <- params[ivpnames]
      names(x) <- sub("\\.0$","",ivpnames)
      x
    }
  }

  if (missing(PACKAGE)) PACKAGE <- character(0)
  
  if (!is.function(rprocess))
    stop("pomp error: 'rprocess' must be a function")
  if (!is.function(dprocess))
    stop("pomp error: 'dprocess' must be a function")

  if (is.function(rmeasure)) {
    if (!all(c('x','t','params','...')%in%names(formals(rmeasure))))
      stop("'rmeasure' must be a function of prototype 'rmeasure(x,t,params,...)'")
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
    stop("'rmeasure' must be either a function or the name of a compiled routine")
  }
  
  if (is.function(dmeasure)) {
    if (!all(c('y','x','t','params','log','...')%in%names(formals(dmeasure))))
      stop("'dmeasure' must be a function of prototype 'dmeasure(y,x,t,params,log,...)'")
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
    stop("'dmeasure' must be either a function or the name of a compiled routine")
  }
  
  if (is.function(skeleton)) {
    if (!all(c('x','t','params','...')%in%names(formals(skeleton))))
      stop("'skeleton' must be a function of prototype 'skeleton(x,t,params,...)'")
    skeleton <- new(
                    "pomp.fun",
                    R.fun=skeleton,
                    use=as.integer(1)
                    )
  } else if (is.character(skeleton)) {
    skeleton <- new(
                    "pomp.fun",
                    R.fun=function(x,t,params,...)stop("very bad: skel.fun"),
                    native.fun=skeleton,
                    PACKAGE=PACKAGE,
                    use=as.integer(2)
                    )
  } else {
    stop("'skeleton' must be either a function or the name of a compiled routine")
  }
  
  if (!is.function(initializer))
    stop("pomp error: 'initializer' must be a function")
  
  if (!all(c('xstart','times','params','...')%in%names(formals(rprocess))))
    stop("pomp error: 'rprocess' must be a function of prototype 'rprocess(xstart,times,params,...)'")
  if (!all(c('x','times','params','log','...')%in%names(formals(dprocess))))
    stop("pomp error: 'dprocess' must be a function of prototype 'dprocess(x,times,params,log,...)'")
  if (!all(c('params','t0','...')%in%names(formals(initializer))))
    stop("pomp error: 'initializer' must be a function of prototype 'initializer(params,t0,...)'")
  
  if (missing(statenames)) statenames <- character(0)
  if (missing(paramnames)) paramnames <- character(0)
  if (missing(covarnames)) covarnames <- character(0)
  
  if (missing(covar)) {
    covar <- matrix(data=0,nrow=0,ncol=0)
    tcovar <- numeric(0)
    covarnames <- character(0)
  } else if (missing(tcovar)) {
    stop("if 'covar' is supplied, 'tcovar' must be supplied also")
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
      stop("if 'covar' is a data frame, 'tcovar' should indicate the time variable")
    } else {
      tcovar <- covar[[tcovar]]
      covar <- as.matrix(covar[-tcovar])
    }
  } else {
    covar <- as.matrix(covar)
  }
  
  if (length(tcovar)!=nrow(covar)) {
    stop("the length of 'tcovar' should match the number of rows of 'covar'")
  } else if (!all(covarnames%in%colnames(covar))) {
    missing <- covarnames[!(covarnames%in%colnames(covar))]
    stop("covariate(s) ",paste(missing,collapse=",")," are not found among the columns of 'covar'")
  } else if (!is.numeric(tcovar)) {
    stop("'tcovar' must either be a numeric vector or must name a numeric vector in the data frame 'covar'")
  }

  if (nrow(covar)>0) {
    if (
        (skeleton@use==1)
        &&!("covars"%in%names(formals(skeleton@R.fun)))
        )
      warning("a covariate table has been given, yet the 'skeleton' function does not have 'covars' as a formal argument")
    if (
        (rmeasure@use==1)
        &&!("covars"%in%names(formals(rmeasure@R.fun)))
        )
      warning("a covariate table has been given, yet the 'rmeasure' function does not have 'covars' as a formal argument")
    if (
        (dmeasure@use==1)
        &&!("covars"%in%names(formals(dmeasure@R.fun)))
        )
      warning("a covariate table has been given, yet the 'dmeasure' function does not have 'covars' as a formal argument")
  }

  new(
      'pomp',
      rprocess = rprocess,
      dprocess = dprocess,
      dmeasure = dmeasure,
      rmeasure = rmeasure,
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
    stop("you must provide at least one formula")
  for (k in 1:nobs) {
    if (!inherits(formulae[[k]],"formula"))
      stop("'measurement.model' takes formulae as arguments")
  }
  obsnames <- unlist(lapply(formulae,function(x)x[[2]]))
  distrib <- lapply(formulae,function(x)as.character(x[[3]][[1]]))
  ddistrib <- lapply(distrib,function(x)paste("d",x,sep=''))
  rdistrib <- lapply(distrib,function(x)paste("r",x,sep=''))
  for (k in 1:nobs) {
    res <- try(
               match.fun(ddistrib[[k]]),
               silent=TRUE
               )
    if (inherits(res,'try-error'))
      stop("distribution function ",ddistrib[[k]]," not found")
    res <- try(
               match.fun(rdistrib[[k]]),
               silent=TRUE
               )
    if (inherits(res,'try-error'))
      stop("random deviate function ",rdistrib[[k]]," not found")
  }
  pred.args <- lapply(formulae,function(x)as.list(x[[3]][-1]))
  dcalls <- vector(mode='list',length=nobs)
  rcalls <- vector(mode='list',length=nobs)
  for (k in 1:nobs) {
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
         for (k in 1:nobs) {
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
         for (k in 1:nobs) {
           y[k] <- eval(
                        rcalls[[k]],
                        envir=as.list(c(x,params,covars,t=t))
                        )
         }
         y
       }
       )
}

