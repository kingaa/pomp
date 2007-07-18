## constructor of the pomp class
pomp <- function (data, times, t0, rprocess, dprocess, rmeasure, dmeasure, ...) {
  if (!is.array(data) || !is.numeric(data))
    stop("'data' must be specified as a numeric array")
  if (!is.numeric(times) || !all(diff(times)>0))
    stop("'times' must be an increasing numeric vector")
  if (length(times)!=ncol(data))
    stop("the length of 'times' does not equal the number of columns in 'data'")
  if (!is.numeric(t0) || length(t0) > 1)
    stop("the zero-time 't0' must be a single number")
  if (t0 >= times[1])
    stop("the zero-time 't0' must occur before the first observation")
  if (missing(rprocess))
    rprocess <- function(xstart,times,params,...)stop("'rprocess' not specified")
  if (missing(dprocess))
    dprocess <- function(x,times,params,log=FALSE,...)stop("'dprocess' not specified")
  if (missing(rmeasure))
    rmeasure <- function(x,times,params,...)stop("'rmeasure' not specified")
  if (missing(dmeasure))
    dmeasure <- function(y,x,times,params,log=FALSE,...)stop("'dmeasure' not specified")
  if (!is.function(rprocess))
    stop("'rprocess' must be a function")
  if (!is.function(dprocess))
    stop("'dprocess' must be a function")
  if (!is.function(rmeasure))
    stop("'rmeasure' must be a function")
  if (!is.function(dmeasure))
    stop("'dmeasure' must be a function")
  if (!all(c('xstart','times','params','...')%in%names(formals(rprocess))))
    stop("'rprocess' must be a function of prototype 'rprocess(xstart,times,params,...)'")
  if (!all(c('x','times','params','log','...')%in%names(formals(dprocess))))
    stop("'dprocess' must be a function of prototype 'dprocess(x,times,params,log,...)'")
  if (!all(c('x','times','params','...')%in%names(formals(rmeasure))))
    stop("'rmeasure' must be a function of prototype 'rmeasure(x,times,params,...)'")
  if (!all(c('y','x','times','params','log','...')%in%names(formals(dmeasure))))
    stop("'dmeasure' must be a function of prototype 'dmeasure(y,x,times,params,log,...)'")
  new(
      'pomp',
      rprocess = rprocess,
      dprocess = dprocess,
      rmeasure = rmeasure,
      dmeasure = dmeasure,
      data = data,
      times = times,
      t0 = t0,
      userdata = list(...)
      )
}
