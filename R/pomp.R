## constructor of the pomp class
pomp <- function (data, times, t0, rprocess, dprocess, rmeasure, dmeasure, initializer, ...) {
  if (is.data.frame(data)) {
    if (!is.character(times) || length(times)!=1 || !(times%in%names(data)))
      stop("pomp error: 'times' must be the name of a column of 'data'")
    tmnm <- times
    times <- data[[tmnm]]
    data <- do.call(rbind,lapply(data[!(names(data)%in%tmnm)],as.numeric))
  }
  if (!is.array(data) || !is.numeric(data))
    stop("pomp error: 'data' must be specified as a numeric array")
  if (!is.numeric(times) || !all(diff(times)>0))
    stop("pomp error: 'times' must be an increasing numeric vector")
  if (length(times)!=ncol(data))
    stop("pomp error: the length of 'times' does not equal the number of columns in 'data'")
  if (!is.numeric(t0) || length(t0) > 1)
    stop("pomp error: the zero-time 't0' must be a single number")
  if (t0 > times[1])
    stop("pomp error: the zero-time 't0' must occur no later than the first observation")
  if (missing(rprocess))
    rprocess <- function(xstart,times,params,...)stop("'rprocess' not specified")
  if (missing(dprocess))
    dprocess <- function(x,times,params,log=FALSE,...)stop("'dprocess' not specified")
  if (missing(rmeasure))
    rmeasure <- function(x,times,params,...)stop("'rmeasure' not specified")
  if (missing(dmeasure))
    dmeasure <- function(y,x,times,params,log=FALSE,...)stop("'dmeasure' not specified")
  if (missing(initializer))
    initializer <- default.initializer
  if (!is.function(rprocess))
    stop("pomp error: 'rprocess' must be a function")
  if (!is.function(dprocess))
    stop("pomp error: 'dprocess' must be a function")
  if (!is.function(rmeasure))
    stop("pomp error: 'rmeasure' must be a function")
  if (!is.function(dmeasure))
    stop("pomp error: 'dmeasure' must be a function")
  if (!is.function(initializer))
    stop("pomp error: 'initializer' must be a function")
  if (!all(c('xstart','times','params','...')%in%names(formals(rprocess))))
    stop("pomp error: 'rprocess' must be a function of prototype 'rprocess(xstart,times,params,...)'")
  if (!all(c('x','times','params','log','...')%in%names(formals(dprocess))))
    stop("pomp error: 'dprocess' must be a function of prototype 'dprocess(x,times,params,log,...)'")
  if (!all(c('x','times','params','...')%in%names(formals(rmeasure))))
    stop("pomp error: 'rmeasure' must be a function of prototype 'rmeasure(x,times,params,...)'")
  if (!all(c('y','x','times','params','log','...')%in%names(formals(dmeasure))))
    stop("pomp error: 'dmeasure' must be a function of prototype 'dmeasure(y,x,times,params,log,...)'")
  if (!all(c('params','t0','...')%in%names(formals(initializer))))
    stop("pomp error: 'initializer' must be a function of prototype 'initializer(params,t0,...)'")
  storage.mode(data) <- 'double'
  storage.mode(times) <- 'double'
  storage.mode(t0) <- 'double'
  new(
      'pomp',
      rprocess = rprocess,
      dprocess = dprocess,
      rmeasure = rmeasure,
      dmeasure = dmeasure,
      data = data,
      times = times,
      t0 = t0,
      initializer = initializer,
      userdata = list(...)
      )
}

default.initializer <- function (params,t0,...) {
  ivpnames <- grep("\\.0$",names(params),val=TRUE)
  if (length(ivpnames)<1)
    stop("no initial value parameters (names ending in '.0') found: see 'pomp' documentation")
  x <- params[ivpnames]
  names(x) <- gsub("\\.0$","",ivpnames)
  x
}
