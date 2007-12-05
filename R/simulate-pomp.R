## simulate a partially-observed Markov process
simulate.pomp.default <- function (object, nsim = 1, seed = NULL, xstart, params,
                                   states = FALSE, obs = FALSE,
                                   times = c(object@t0,time(object)), ...) {
  ntimes <- length(times)
  times <- as.numeric(times)
  if (ntimes<1)
    stop("if length of 'times' is less than 1, there is no work to do")
  if (is.null(dim(xstart)))
    xstart <- matrix(xstart,ncol=1,dimnames=list(names(xstart),NULL))
  if (is.null(dim(params)))
    params <- matrix(params,ncol=1,dimnames=list(names(params),NULL))
  if (is.null(rownames(xstart)))
    stop("'xstart' must have rownames")
  if (is.null(rownames(params)))
    stop("'params' must have rownames")
  if ((ncol(xstart)==1)&&(ncol(params)>1))
    xstart <- matrix(xstart,nrow=nrow(xstart),ncol=ncol(params),dimnames=list(rownames(xstart),NULL))
  if ((ncol(params)==1)&&(ncol(xstart)>1))
    params <- matrix(params,nrow=nrow(params),ncol=ncol(xstart),dimnames=list(rownames(params),NULL))
  npars <- ncol(xstart)
  if (npars!=ncol(params))
    stop("'xstart' and 'params' must have equal number of columns")
  if (!is.null(seed)) { # set the random seed (be very careful about this)
    if (!exists('.Random.seed',envir=.GlobalEnv)) runif(1)
    save.seed <- get('.Random.seed',envir=.GlobalEnv)
    set.seed(seed)
  }
  nreps <- npars*nsim                   # total number of replicates
  ## we will do the process model simulations with single calls to the user functions
  if (nsim > 1) {    # make nsim copies of the IC and parameter matrices
    xstart <- array(rep(xstart,nsim),dim=c(nrow(xstart),nreps),dimnames=list(rownames(xstart),NULL))
    params <- array(rep(params,nsim),dim=c(nrow(params),nreps),dimnames=list(rownames(params),NULL))
  }
  x <- rprocess(object,xstart=xstart,times=times,params=params) # simulate the process model
  ## rprocess returns a rank-3 matrix (nvars x nreps x ntimes)  
  if (obs || !states) { # we need only simulate the observation process if obs=T or states=F
    y <- rmeasure(object,x=x,times=times,params=params)
    ## rmeasure returns a rank-3 matrix (nobs x nreps x ntimes)
  }
  if (!is.null(seed)) {                 # restore the RNG state
    assign('.Random.seed',save.seed,envir=.GlobalEnv)
  }
  if (!obs && !states) { # both obs=F and states=F, return a list of pomp objects
    nm.y <- rownames(y)
    nobs <- nrow(y)
    retval <- lapply(
                     seq(length=nreps),
                     function (rep) {
                       x <- as(object,'pomp') # copy the pomp object
                       x@t0 <- times[1]
                       x@times <- times[-1]
                       x@data <- array(0,dim=c(nobs,ntimes-1),dimnames=list(nm.y,NULL))
                       x@data[,] <- y[,rep,-1] # replace the data
                       x
                     }
                     )
  } else {
    if (states) {                       # reformat the x array
      nvars <- nrow(x)
      nm.x <- rownames(x)
      dim(x) <- c(nvars,npars,nsim,ntimes)
      rownames(x) <- nm.x
    }
    if (obs) {                          # reformat the y array
      nobs <- nrow(y)
      nm.y <- rownames(y)
      dim(y) <- c(nobs,npars,nsim,ntimes)
      rownames(y) <- nm.y
    }
    if (!obs && states) {               # return states only
      retval <- x
    } else if (obs && !states) {        # return only observations
      retval <- y
    } else {                     # return both states and observations
      retval <- list(
                     states = x,
                     obs = y
                     )
    }
  }
  retval
}

setMethod('simulate','pomp',simulate.pomp.default)
