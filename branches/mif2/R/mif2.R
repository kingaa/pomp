## MIF algorithm functions

## define the mif class
setClass(
         'mif2d.pomp',
         contains='pfilterd.pomp',
         representation=representation(
           transform = "logical",
           Nmif = 'integer',
           perturb.fn = 'function',
           conv.rec = 'matrix'
           )
         )


mif2.pfilter <- function (object, params, Np,
                          tol, max.fail,
                          pred.mean, pred.var, filter.mean,
                          mifiter, perturb.fn,
                          transform, verbose,
                          .getnativesymbolinfo = TRUE) {
  
  ptsi.inv <- ptsi.for <- gnsi.rproc <- gnsi.dmeas <- as.logical(.getnativesymbolinfo)
  transform <- as.logical(transform)
  mifiter <- as.integer(mifiter)
  Np <- as.integer(Np)
  Np <- c(Np,Np[1])
  
  times <- time(object,t0=TRUE)
  ntimes <- length(times)-1
  
  paramnames <- rownames(params)
  npars <- nrow(params)

  ## perturb parameters
  params <- apply(params,2L,perturb.fn,mifiter=mifiter,timeno=1L)
    
  ## transform parameters if necessary
  if (transform) {
    tparams <- partrans(object,params,dir="forward",
                        .getnativesymbolinfo=ptsi.for)
    ptsi.for <- FALSE
  }

  ## get initial states
  x <- init.state(
                  object,
                  params=if (transform) tparams else params,
                  )
  statenames <- rownames(x)
  nvars <- nrow(x)
  
  loglik <- rep(NA,ntimes)
  eff.sample.size <- numeric(ntimes)
  nfail <- 0
  
  ## set up storage for prediction means, variances, etc.
  if (pred.mean)
    pred.m <- matrix(
                     data=0,
                     nrow=nvars+npars,
                     ncol=ntimes,
                     dimnames=list(c(statenames,paramnames),NULL)
                     )
  else
    pred.m <- array(data=numeric(0),dim=c(0,0))
  
  if (pred.var)
    pred.v <- matrix(
                     data=0,
                     nrow=nvars+npars,
                     ncol=ntimes,
                     dimnames=list(c(statenames,paramnames),NULL)
                     )
  else
    pred.v <- array(data=numeric(0),dim=c(0,0))
  
  if (filter.mean)
    filt.m <- matrix(
                     data=0,
                     nrow=nvars+length(paramnames),
                     ncol=ntimes,
                     dimnames=list(c(statenames,paramnames),NULL)
                     )
  else
    filt.m <- array(data=numeric(0),dim=c(0,0))

  for (nt in seq_len(ntimes)) {

    if (nt > 1) {
      ## perturb parameters
      params <- apply(params,2L,perturb.fn,mifiter=mifiter,timeno=nt)
      ## transform parameters if necessary
      if (transform) {
        tparams <- partrans(object,params,dir="forward",
                            .getnativesymbolinfo=ptsi.for)
      }
    }
    
    ## advance the state variables according to the process model
    X <- try(
             rprocess(
                      object,
                      xstart=x,
                      times=times[c(nt,nt+1)],
                      params=if (transform) tparams else params,
                      offset=1,
                      .getnativesymbolinfo=gnsi.rproc
                      ),
             silent=FALSE
             )
    if (inherits(X,'try-error'))
      stop(sQuote("mif2.pfilter")," error: process simulation error")
    gnsi.rproc <- FALSE
    
    ## determine the weights
    weights <- try(
                   dmeasure(
                            object,
                            y=object@data[,nt,drop=FALSE],
                            x=X,
                            times=times[nt+1],
                            params=if (transform) tparams else params,
                            log=FALSE,
                            .getnativesymbolinfo=gnsi.dmeas
                            ),
                   silent=FALSE
                   )
    if (inherits(weights,'try-error'))
      stop(sQuote("mif2.pfilter")," error: error in calculation of weights")
    if (any(!is.finite(weights))) {
      stop(sQuote("mif2.pfilter")," error: ",sQuote("dmeasure"),
           " returns non-finite value")
    }
    gnsi.dmeas <- FALSE
    
    ## compute prediction mean, prediction variance, filtering mean,
    ## effective sample size, log-likelihood
    ## also do resampling if filtering has not failed
    xx <- try(
              .Call(
                    pomp:::pfilter_computations,
                    X,params,Np[nt+1],
                    FALSE,numeric(0),
                    pred.mean,pred.var,filter.mean,
                    FALSE,weights,tol
                    ),
              silent=FALSE
              )
    if (inherits(xx,'try-error')) {
      stop(sQuote("mif2.pfilter")," error",call.=FALSE)
    }
    all.fail <- xx$fail
    loglik[nt] <- xx$loglik
    eff.sample.size[nt] <- xx$ess
    
    x <- xx$states
    params <- xx$params
    
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
      if (nfail>max.fail)
        stop(sQuote("mif2.pfilter")," error: too many filtering failures",call.=FALSE)
    }
    
    if (verbose && (nt%%5==0))
      cat("mif2 pfilter timestep",nt,"of",ntimes,"finished\n")
    
  }
  
  if (nfail>0)
    warning(sprintf(ngettext(nfail,msg1="%d filtering failure occurred in ",
                             msg2="%d filtering failures occurred in "),nfail),
            sQuote("mif2.pfilter"),call.=FALSE)

  new(
      "pfilterd.pomp",
      object,
      pred.mean=pred.m,
      pred.var=pred.v,
      filter.mean=filt.m,
      paramMatrix=params,
      eff.sample.size=eff.sample.size,
      cond.loglik=loglik,
      saved.states=list(),
      saved.params=list(),
      seed=as.integer(NULL),
      Np=as.integer(head(Np,-1)),
      tol=tol,
      nfail=as.integer(nfail),
      loglik=sum(loglik)
      )
}

mif2.internal <- function (object, Nmif,
                           start, Np, perturb.fn,
                           tol, max.fail,
                           verbose, transform, .ndone = 0L,
                           .getnativesymbolinfo = TRUE,
                           ...) {
  
  gnsi <- as.logical(.getnativesymbolinfo)
  transform <- as.logical(transform)
  
  Nmif <- as.integer(Nmif)
  if (Nmif<0)
    stop("mif error: ",sQuote("Nmif")," must be a positive integer",call.=FALSE)

  if (transform)
    start <- partrans(object,start,dir="inverse")
  
  ntimes <- length(time(object))

  start.names <- names(start)
  if (is.null(start.names))
    stop("mif2 error: ",sQuote("start")," must be a named vector")
  
  if (!is.function(perturb.fn))
    stop("mif2 error: ",sQuote("perturb.fn")," must be a function")

  if (is.function(Np)) {
    Np <- try(
              vapply(seq.int(from=1,to=ntimes,by=1),Np,numeric(1)),
              silent=FALSE
              )
    if (inherits(Np,"try-error"))
      stop("if ",sQuote("Np")," is a function, it must return a single positive integer")
  }
  if (length(Np)==1)
    Np <- rep(Np,times=ntimes)
  else if (length(Np)!=ntimes)
    stop(sQuote("Np")," must have length 1 or length ",ntimes)
  if (any(Np<=0))
    stop("number of particles, ",sQuote("Np"),", must always be positive")
  if (!is.numeric(Np))
    stop(sQuote("Np")," must be a number, a vector of numbers, or a function")
  Np <- as.integer(Np)
  
  conv.rec <- matrix(
                     data=NA,
                     nrow=Nmif+1,
                     ncol=length(start)+2,
                     dimnames=list(
                       seq(.ndone,.ndone+Nmif),
                       c('loglik','nfail',names(start))
                       )
                     )
  conv.rec[1L,] <- c(NA,NA,start)

  if (.ndone > 0) {                     # call is from 'continue'
    paramMatrix <- object@paramMatrix
  } else if (Nmif > 0) {                # initial call
    paramMatrix <- array(data=start,dim=c(length(start),Np[1L]),
                         dimnames=list(names(start),NULL))
  } else {                              # no work to do
    paramMatrix <- array(data=numeric(0),dim=c(0,0))
  }

  object <- as(object,"pomp")
    
  for (n in seq_len(Nmif)) { ## iterate the filtering

    pfp <- try(
               mif2.pfilter(
                            object=object,
                            params=paramMatrix,
                            Np=Np,
                            tol=tol,
                            max.fail=max.fail,
                            pred.mean=(n==Nmif),
                            pred.var=(n==Nmif),
                            filter.mean=(n==Nmif),
                            mifiter=.ndone+n,
                            perturb.fn=perturb.fn,
                            transform=transform,
                            verbose=verbose,
                            .getnativesymbolinfo=gnsi
                            ),
               silent=FALSE
               )
    if (inherits(pfp,"try-error")) 
      stop("mif2 particle-filter error")

    gnsi <- FALSE
    
    theta <- rowMeans(pfp@paramMatrix[,,drop=FALSE])

    conv.rec[n+1,-c(1,2)] <- theta
    conv.rec[n,c(1,2)] <- c(pfp@loglik,pfp@nfail)
    
    if (verbose) cat("MIF2 iteration ",n," of ",Nmif," completed\n")
    
  } ### end of main loop

  ## back transform the parameter estimate if necessary
  if (transform) theta <- partrans(pfp,theta,dir="forward")
  
  new(
      "mif2d.pomp",
      pfp,
      params=theta,
      paramMatrix=pfp@paramMatrix,
      transform=transform,
      Nmif=Nmif,
      perturb.fn=perturb.fn,
      conv.rec=conv.rec,
      tol=tol
      )
}

setMethod(
          "mif2",
          signature=signature(object="pomp"),
          function (object, Nmif = 1,
                    start, Np, perturb.fn,
                    tol = 1e-17, max.fail = Inf,
                    verbose = getOption("verbose"),
                    transform = FALSE,
                    ...) {
            
            if (missing(start)) start <- coef(object)
            if (length(start)==0)
              stop(
                   "mif2 error: ",sQuote("start")," must be specified if ",
                   sQuote("coef(object)")," is NULL",
                   call.=FALSE
                   )

            if (missing(Np))
              stop("mif2 error: ",sQuote("Np")," must be specified",call.=FALSE)

            if (missing(perturb.fn))
              stop("mif2 error: ",sQuote("perturb.fn")," must be specified",call.=FALSE)
            perturb.fn <- match.fun(perturb.fn)
            if (!all(c('params','mifiter','timeno','...')%in%names(formals(perturb.fn)))) {
              stop(
                   "mif2 error: ",
                   sQuote("perturb.fn"),
                   " must be a function of prototype ",
                   sQuote("perturb.fn(params,mifiter,timeno,...)"),
                   call.=FALSE
                   )
            }
            
            mif2.internal(
                          object=object,
                          Nmif=Nmif,
                          start=start,
                          Np=Np,
                          perturb.fn=perturb.fn,
                          tol=tol,
                          max.fail=max.fail,
                          transform=transform,
                          verbose=verbose,
                          ...
                          )
            
          }
          )


setMethod(
          "mif2",
          signature=signature(object="pfilterd.pomp"),
          function (object, Nmif = 1, Np, tol, ...) {
            
            if (missing(Np)) Np <- object@Np
            if (missing(tol)) tol <- object@tol

            mif2(
                 object=as(object,"pomp"),
                 Nmif=Nmif,
                 Np=Np,
                 tol=tol,
                 ...
                 )
          }
          )

setMethod(
          "mif2",
          signature=signature(object="mif2d.pomp"),
          function (object, Nmif, start,
                    Np, perturb.fn,
                    tol, transform,
                    ...) {
            
            if (missing(Nmif)) Nmif <- object@Nmif
            if (missing(start)) start <- coef(object)
            if (missing(perturb.fn)) perturb.fn <- object@perturb.fn
            if (missing(transform)) transform <- object@transform

            if (missing(Np)) Np <- object@Np
            if (missing(tol)) tol <- object@tol

            mif2(
                 object=as(object,"pomp"),
                 Nmif=Nmif,
                 start=start,
                 Np=Np,
                 perturb.fn=perturb.fn,
                 tol=tol,
                 transform=transform,
                 ...
                 )
          }
          )

setMethod(
          'continue',
          signature=signature(object='mif2d.pomp'),
          function (object, Nmif = 1,
                    ...) {
            
            ndone <- object@Nmif
            
            obj <- mif2(
                        object=object,
                        Nmif=Nmif,
                        .ndone=ndone,
                        ...
                        )
            
            object@conv.rec[ndone+1,c('loglik','nfail')] <- obj@conv.rec[1L,c('loglik','nfail')]
            obj@conv.rec <- rbind(
                                  object@conv.rec,
                                  obj@conv.rec[-1L,colnames(object@conv.rec)]
                                  )
            obj@Nmif <- as.integer(ndone+Nmif)
            
            obj
          }
          )

default.cooling <- function (object, fraction, par.sd, ic.sd) {
  
  nT <- length(time(object))
  theta <- (1-fraction)/fraction/(50*nT-1)

  function (params, mifiter, timeno, ...) {
    pert <- params
    sigma <- 1/(1+theta*((mifiter-1)*nT+timeno-1))
    if (timeno==1) {
      pert[names(ic.sd)] <- rnorm(
                                  n=length(ic.sd),
                                  mean=pert[names(ic.sd)],
                                  sd=ic.sd*sigma
                                  )
    }
    pert[names(par.sd)] <- rnorm(
                                 n=length(par.sd),
                                 mean=pert[names(par.sd)],
                                 sd=par.sd*sigma
                                 )
    pert
  }
}
