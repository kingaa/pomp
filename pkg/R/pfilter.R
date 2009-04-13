## particle filtering codes

## question: when pfilter.internal is called by mif, do we need to compute the prediction means and variances of the state variables each time, or only at the end?
## question: how much efficiency would be realized by eliminating the calls to 'apply' with something else?

pfilter.internal <- function (object, params, Np,
                              tol, warn, max.fail,
                              pred.mean, pred.var, filter.mean,
                              .rw.sd, verbose) {
  if (missing(params)) {
    params <- coef(object)
    if (length(params)==0) {
      stop("pfilter error: ",sQuote("params")," must be supplied",call.=FALSE)
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
    stop("pfilter error: ",sQuote("params")," must have rownames",call.=FALSE)

  x <- init.state(object,params=params)
  statenames <- rownames(x)
  nvars <- nrow(x)
  
  random.walk <- !missing(.rw.sd)
  if (random.walk) {
    rw.names <- names(.rw.sd)
    if (is.null(rw.names)||!is.numeric(.rw.sd))
      stop("pfilter error: ",sQuote(".rw.sd")," must be a named vector",call.=FALSE)
    if (any(!(rw.names%in%paramnames)))
      stop("pfilter error: the rownames of ",sQuote("params")," must include all of the names of ",sQuote(".rw.sd"),"",call.=FALSE)
    sigma <- .rw.sd
  } else {
    rw.names <- character(0)
  }
  
  loglik <- rep(NA,ntimes)
  eff.sample.size <- rep(NA,ntimes)
  nfail <- 0
  npars <- length(rw.names)
  
  pred.m <-  if (pred.mean)
    matrix(
           data=0,
           nrow=nvars+npars,
           ncol=ntimes,
           dimnames=list(c(statenames,rw.names),NULL)
           )
  else NULL
  
  pred.v <- if (pred.var)
    matrix(
           data=0,
           nrow=nvars+npars,
           ncol=ntimes,
           dimnames=list(c(statenames,rw.names),NULL)
           )
  else NULL
  
  filt.m <- if (filter.mean)
    if (random.walk) 
      matrix(
             data=0,
             nrow=nvars+length(paramnames),
             ncol=ntimes,
             dimnames=list(c(statenames,paramnames),NULL)
             )
    else
      matrix(
             data=0,
             nrow=nvars,
             ncol=ntimes,
             dimnames=list(statenames,NULL)
             )
  else NULL

  for (nt in seq(length=ntimes)) {
    
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
      stop("pfilter error: process simulation error",call.=FALSE)

    x[,] <- X                 # ditch the third dimension
    
    ## prediction means
    if (pred.mean) {                    
      xx <- try(
                c(
                  rowMeans(x),
                  rowMeans(params[rw.names,,drop=FALSE])
                  ),
                silent=FALSE
                )
      if (inherits(xx,'try-error')) {
        stop("pfilter error: error in prediction mean computation",call.=FALSE)
      } else {
        pred.m[,nt] <- xx
      }
    }

    ## prediction variances
    if (pred.var) {
      problem.indices <- unique(which(!is.finite(x),arr.ind=TRUE)[,1])
      if (length(problem.indices)>0) {
        stop(
             "pfilter error: non-finite state variable(s): ",
             paste(rownames(x)[problem.indices],collapse=', '),
             call.=FALSE
             )
      }
      if (random.walk) {
        problem.indices <- unique(which(!is.finite(params[rw.names,,drop=FALSE]),arr.ind=TRUE)[,1])
        if (length(problem.indices)>0) {
          stop(
               "pfilter error: non-finite parameter(s): ",
               paste(rw.names[problem.indices],collapse=', '),
               call.=FALSE
               )
        }
      }
      xx <- try(
                c(
                  apply(x,1,var),
                  apply(params[rw.names,,drop=FALSE],1,var)
                  ),
                silent=FALSE
                )
      if (inherits(xx,'try-error')) {
        stop("pfilter error: error in prediction variance computation",call.=FALSE)
      } else {
        pred.v[,nt] <- xx
      }
    }

    ## determine the weights
    weights <- try(
                   dmeasure(
                            object,
                            y=object@data[,nt,drop=FALSE],
                            x=X,
                            times=times[nt+1],
                            params=params
                            ),
                   silent=FALSE
                   )
    if (inherits(weights,'try-error'))
      stop("pfilter error: error in calculation of weights",call.=FALSE)
    if (any(is.na(weights))) {
      ## problem.indices <- which(is.na(weights))
      stop("pfilter error: dmeasure returns NA",call.=FALSE)
    }

    ## test for failure to filter
    dim(weights) <- NULL
    failures <- weights < tol
    all.fail <- all(failures)
    if (all.fail) {                     # all particles are lost
      if (warn)
        message("filtering failure at time t = ",times[nt+1])
      nfail <- nfail+1
      if (nfail > max.fail)
        stop('pfilter error: too many filtering failures',call.=FALSE)
      loglik[nt] <- log(tol)          # worst log-likelihood
      weights <- rep(1/Np,Np)
      eff.sample.size[nt] <- 0
    } else {                  # not all particles are lost
      ## compute log-likelihood
      loglik[nt] <- log(mean(weights))  
      weights[failures] <- 0
      weights <- weights/sum(weights)
      ## compute effective sample-size
      eff.sample.size[nt] <- 1/(weights%*%weights) 
    }

    ## compute filtering means
    if (filter.mean) {
      filt.m[statenames,nt] <- x %*% weights
      if (random.walk)
        filt.m[paramnames,nt] <- params %*% weights
    }

    ## Matrix with samples (columns) from filtering distribution theta.t | Y.t
    if (!all.fail) {
      sample <- .Call(systematic_resampling,weights)
      x <- x[,sample,drop=FALSE]
      params <- params[,sample,drop=FALSE]
    }

    ## random walk for parameters
    if (random.walk) {
      pred.v[rw.names,nt] <- pred.v[rw.names,nt]+sigma^2
      params[rw.names,] <- params[rw.names,]+rnorm(n=Np*length(sigma),mean=0,sd=sigma)
    }

    if (verbose && ((ntimes-nt)%%5==0))
      cat("pfilter timestep",nt,"of",ntimes,"finished\n")

  }

  list(
       pred.mean=pred.m,
       pred.var=pred.v,
       filter.mean=filt.m,
       eff.sample.size=eff.sample.size,
       cond.loglik=loglik,
       nfail=nfail,
       loglik=sum(loglik)
       )
}

setMethod(
          "pfilter",
          "pomp",
          function (object, params, Np,
                    tol = 1e-17, warn = TRUE, max.fail = 0,
                    pred.mean = FALSE,
                    pred.var = FALSE,
                    filter.mean = FALSE,
                    verbose = FALSE, ...) {
            pfilter.internal(
                             object=object,
                             params=params,
                             Np=Np,
                             tol=tol,
                             warn=warn,
                             max.fail=max.fail,
                             pred.mean=pred.mean,
                             pred.var=pred.var,
                             filter.mean=filter.mean,
                             verbose=verbose
                             )
          }
          )
