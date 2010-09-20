## particle filtering codes

## generic particle filter
pfilter <- function (object, ...)
  stop("function ",sQuote("pfilter")," is undefined for objects of class ",sQuote(class(object)))
setGeneric('pfilter')  

## question: when pfilter.internal is called by mif, do we need to compute the prediction means and variances of the state variables each time, or only at the end?
## question: how much efficiency would be realized by eliminating the calls to 'apply' with something else?

pfilter.internal <- function (object, params, Np,
                              tol, max.fail,
                              pred.mean, pred.var, filter.mean,
                              .rw.sd, seed, verbose,
                              save.states) {
  if (missing(seed)) seed <- NULL
  if (!is.null(seed)) {
    if (!exists(".Random.seed",where=.GlobalEnv)) { # need to initialize the RNG
      runif(1)
    }
    save.seed <- get(".Random.seed",pos=.GlobalEnv)
    set.seed(seed)
  }

  if (missing(params)) {
    params <- coef(object)
    if (length(params)==0) {
      stop(sQuote("pfilter")," error: ",sQuote("params")," must be supplied",call.=FALSE)
    }
  }
  if (missing(Np))
    Np <- NCOL(params)
  times <- time(object,t0=TRUE)
  ntimes <- length(times)-1
  if (is.null(dim(params))) {
    params <- matrix(
                     params,
                     nrow=length(params),
                     ncol=Np,
                     dimnames=list(
                       names(params),
                       NULL
                       )
                     )
  }
  paramnames <- rownames(params)
  if (is.null(paramnames))
    stop(sQuote("pfilter")," error: ",sQuote("params")," must have rownames",call.=FALSE)

  x <- init.state(object,params=params)
  statenames <- rownames(x)
  nvars <- nrow(x)
  if (save.states) {
    xparticles <- array(
                        data=NA,
                        dim=c(nvars,Np,ntimes),
                        dimnames=list(statenames,NULL,NULL)
                        )
  } else {
    xparticles <- NULL
  }
  
  random.walk <- !missing(.rw.sd)
  if (random.walk) {
    rw.names <- names(.rw.sd)
    if (is.null(rw.names)||!is.numeric(.rw.sd))
      stop(sQuote("pfilter")," error: ",sQuote(".rw.sd")," must be a named vector",call.=FALSE)
    if (any(!(rw.names%in%paramnames)))
      stop(
           sQuote("pfilter")," error: the rownames of ",
           sQuote("params")," must include all of the names of ",
           sQuote(".rw.sd"),"",call.=FALSE
           )
    sigma <- .rw.sd
  } else {
    rw.names <- character(0)
  }
  
  loglik <- rep(NA,ntimes)
  eff.sample.size <- numeric(ntimes)
  nfail <- 0
  npars <- length(rw.names)

  pred.m <- NULL
  pred.v <- NULL
  filt.m <- NULL

  ## set up storage for prediction means, variances, etc.
  if (pred.mean)
    pred.m <- matrix(
                     data=0,
                     nrow=nvars+npars,
                     ncol=ntimes,
                     dimnames=list(c(statenames,rw.names),NULL)
                     )
  
  if (pred.var)
    pred.v <- matrix(
                     data=0,
                     nrow=nvars+npars,
                     ncol=ntimes,
                     dimnames=list(c(statenames,rw.names),NULL)
                     )
  
  if (filter.mean) {
    if (random.walk) {
      filt.m <- matrix(
                       data=0,
                       nrow=nvars+length(paramnames),
                       ncol=ntimes,
                       dimnames=list(c(statenames,paramnames),NULL)
                       )
    } else {
      filt.m <- matrix(
                       data=0,
                       nrow=nvars,
                       ncol=ntimes,
                       dimnames=list(statenames,NULL)
                       )
    }
  }

  for (nt in seq_len(ntimes)) {
    
    ## advance the state variables according to the process model
    X <- try(
             rprocess(
                      object,
                      x=x,
                      times=times[c(nt,nt+1)],
                      params=params
                      )[,,2,drop=FALSE],
             silent=FALSE
             )
    if (inherits(X,'try-error'))
      stop(sQuote("pfilter")," error: process simulation error",call.=FALSE)

    x[,] <- X                 # ditch the third dimension
    
    if (pred.var) { ## check for nonfinite state variables and parameters
      problem.indices <- unique(which(!is.finite(x),arr.ind=TRUE)[,1])
      if (length(problem.indices)>0) {  # state variables
        stop(
             sQuote("pfilter")," error: non-finite state variable(s): ",
             paste(rownames(x)[problem.indices],collapse=', '),
             call.=FALSE
             )
      }
      if (random.walk) { # parameters (need to be checked only if 'random.walk=TRUE')
        problem.indices <- unique(which(!is.finite(params[rw.names,,drop=FALSE]),arr.ind=TRUE)[,1])
        if (length(problem.indices)>0) {
          stop(
               sQuote("pfilter")," error: non-finite parameter(s): ",
               paste(rw.names[problem.indices],collapse=', '),
               call.=FALSE
               )
        }
      }
    }
    
    ## determine the weights
    weights <- try(
                   dmeasure(
                            object,
                            y=object@data[,nt,drop=FALSE],
                            x=X,
                            times=times[nt+1],
                            params=params,
                            log=FALSE
                            ),
                   silent=FALSE
                   )
    if (inherits(weights,'try-error'))
      stop(sQuote("pfilter")," error: error in calculation of weights",call.=FALSE)
    if (any(is.na(weights))) {
      stop(sQuote("pfilter")," error: ",sQuote("dmeasure")," returns NA",call.=FALSE)
    }

    ## prediction mean, prediction variance, filtering mean, effective sample size, log-likelihood
    xx <- try(
              .Call(
                    pfilter_computations,
                    x,params,
                    random.walk,rw.names,
                    pred.mean,pred.var,
                    filter.mean,weights,tol
                    ),
              silent=FALSE
              )
    if (inherits(xx,'try-error')) {
      stop(sQuote("pfilter")," error: error in prediction mean/variance computation",call.=FALSE)
    }
    all.fail <- xx$fail
    loglik[nt] <- xx$loglik
    eff.sample.size[nt] <- xx$ess

    if (pred.mean)
      pred.m[,nt] <- xx$pm
    if (pred.var)
      pred.v[,nt] <- xx$pv
    if (filter.mean)
      filt.m[,nt] <- xx$fm

    if (all.fail) { ## all particles are lost
      nfail <- nfail+1
      if (verbose)
        message("filtering failure at time t = ",times[nt+1])
      if (nfail > max.fail)
        stop(sQuote("pfilter")," error: too many filtering failures",call.=FALSE)
    } else { ## matrix with samples (columns) from filtering distribution theta.t | Y.t
      sample <- .Call(systematic_resampling,weights)
      x <- x[,sample,drop=FALSE]
      params <- params[,sample,drop=FALSE]
    }

    ## random walk for parameters
    if (random.walk) {
      pred.v[rw.names,nt] <- pred.v[rw.names,nt]+sigma^2
      params[rw.names,] <- rnorm(n=Np*length(sigma),mean=params[rw.names,],sd=sigma)
    }

    if (save.states) {
      xparticles[,,nt] <- x
    }

    if (verbose && ((ntimes-nt)%%5==0))
      cat("pfilter timestep",nt,"of",ntimes,"finished\n")

  }

  if (!is.null(seed)) {
    assign(".Random.seed",save.seed,pos=.GlobalEnv)
    seed <- save.seed
  }

  list(
       pred.mean=pred.m,
       pred.var=pred.v,
       filter.mean=filt.m,
       eff.sample.size=eff.sample.size,
       cond.loglik=loglik,
       states=xparticles,
       seed=seed,
       Np=Np,
       tol=tol,
       nfail=nfail,
       loglik=sum(loglik)
       )
}

setMethod(
          "pfilter",
          "pomp",
          function (object, params, Np,
                    tol = 1e-17,
                    max.fail = 0,
                    pred.mean = FALSE,
                    pred.var = FALSE,
                    filter.mean = FALSE,
                    save.states = FALSE,
                    seed = NULL,
                    verbose = getOption("verbose"),
                    ...) {
            pfilter.internal(
                             object=object,
                             params=params,
                             Np=Np,
                             tol=tol,
                             max.fail=max.fail,
                             pred.mean=pred.mean,
                             pred.var=pred.var,
                             filter.mean=filter.mean,
                             save.states=save.states,
                             seed=seed,
                             verbose=verbose
                             )
          }
          )
