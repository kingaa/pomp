## simulate a partially-observed Markov process

simulate.internal <- function (object, nsim = 1, seed = NULL, params,
                               states = FALSE, obs = FALSE,
                               times, t0, ...) {
  warn.condition <- (missing(t0) && ((obs || states) || (!missing(times))))
  if (warn.condition) 
    warning(
            "The default behavior of ",sQuote("simulate")," has changed.\n",
            "See the documentation (",dQuote("pomp?simulate"),") for details\n",
            "and set the ",sQuote("times")," and ",sQuote("t0"),
            " arguments appropriately to compensate.\n",
            call.=FALSE
            )

  if (missing(times)) {
    times <- time(object,t0=FALSE)
  } else {
    times <- as.numeric(times)
  }

  if (length(times)==0)
    stop("if ",sQuote("times")," is empty, there is no work to do",call.=FALSE)
  
  if (any(diff(times)<0))
    stop(sQuote("times")," must be a nondecreasing sequence of times",call.=FALSE)

  if (missing(t0)) {
    t0 <- timezero(object)
  } else {
    t0 <- as.numeric(t0)
  }
  
  if (t0>times[1])
    stop("the zero-time ",sQuote("t0")," must occur no later than the first observation",call.=FALSE)
  ntimes <- length(times)
  
  if (missing(params)) {
    if (length(object@params)>0)
      params <- object@params
    else
      stop("no ",sQuote("params")," specified",call.=FALSE)
  }
  if (is.null(dim(params)))
    params <- matrix(params,ncol=1,dimnames=list(names(params),NULL))
  npars <- ncol(params)

  if (!is.null(seed)) { # set the random seed (be very careful about this)
    if (!exists('.Random.seed',envir=.GlobalEnv)) runif(1)
    save.seed <- get('.Random.seed',envir=.GlobalEnv)
    set.seed(seed)
  }
  
  if (!is.numeric(nsim)||(length(nsim)>1)||nsim<1)
    stop(sQuote("nsim")," must be a positive integer")
  nsim <- as.integer(nsim)

  xstart <- init.state(object,params=params,t0=t0)
  nreps <- npars*nsim                   # total number of replicates
  ## we will do the process model simulations with single calls to the user functions
  if (nsim > 1) {    # make nsim copies of the IC and parameter matrices
    xstart <- array(rep(xstart,nsim),dim=c(nrow(xstart),nreps),dimnames=list(rownames(xstart),NULL))
    params <- array(rep(params,nsim),dim=c(nrow(params),nreps),dimnames=list(rownames(params),NULL))
  }

  x <- rprocess(object,xstart=xstart,times=c(t0,times),params=params) # simulate the process model
  x <- x[,,-1,drop=FALSE]
  ## rprocess returns a rank-3 matrix (nvars x nreps x (ntimes+1))  

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
                     seq_len(nreps),
                     function (rep) {
                       po <- as(object,'pomp')
                       po@t0 <- t0
                       po@times <- times
                       po@params <- params[,rep]
                       po@data <- array(0,dim=c(nobs,ntimes),dimnames=list(nm.y,NULL))
                       po@data[,] <- y[,rep,]
                       po@states <- array(dim=c(nrow(x),ntimes),dimnames=list(rownames(x),NULL))
                       po@states[,] <- x[,rep,]
                       po
                     }
                     )
    if (nreps==1) retval <- retval[[1]]
  } else if (!obs && states)            # return states only
    retval <- x
  else if (obs && !states)              # return only observations
    retval <- y
  else                           # return both states and observations
    retval <- list(
                   states = x,
                   obs = y
                   )

  retval
}

setMethod("simulate",signature=signature(object="pomp"),definition=simulate.internal)
