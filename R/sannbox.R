######################
## Simulated annealing minimizer with box constraints
##
## By default, the annealing schedule is
## temp / log(((k-1) %/% tmax)*tmax + exp(1)), where
## the parameters of this schedule can be changed via
## the 'control' argument, and k ranges from 0 to 'maxit'
##
## modified from code originally written by 
## Daniel Reuman, Imperial College London

sannbox <- function (par, fn, control = list(), ...) {

  big <- 1e35  ## a very large number

  npar <- length(par)
  neval <- 0
  
  control.default <- list(
                          maxit=10000,
                          temp=1,
                          tmax=10,
                          sched=NULL,
                          candidate.dist=NULL,
                          fnscale=1,
                          parscale=1,
                          lower=-Inf,
                          upper=Inf,
                          trace=0
                          )
  control.default[names(control)] <- control
  control <- control.default

  if (is.null(control$sched))           # default cooling schedule
    control$sched <- function (k, temp, tmax) temp/log(((k-1)%/%tmax)*tmax+exp(1))

  if (is.function(control$sched))
    temps <- vapply(
                    seq_len(control$maxit),
                    FUN=control$sched,
                    FUN.VALUE=numeric(1),
                    temp=control$temp,
                    tmax=control$tmax
                    )
  else if (is.numeric(control$sched)) {
    temps <- control$sched
    if (length(temps)<control$maxit)
      stop("insufficiently many temperatures supplied in ",sQuote("control$sched"))
  }
  
  if (is.null(control$candidate.dist))
    candidate.dist <- function (par, temp, scale) rnorm(n=npar,mean=par,sd=scale*temp)

  if (length(control$lower)<npar)
    control$lower <- rep(control$lower,npar)
  if (length(control$upper)<npar)
    control$upper <- rep(control$upper,npar)

  ## initialization for the algorithm
  thetabest <- thetacurrent <- par
  ycurrent <- fn(thetacurrent,...)/control$fnscale
  if (!is.finite(ycurrent)) ycurrent <- big
  ybest <- ycurrent
  neval <- 1
  
  if (control$trace>0)
    cat("initial evaluation: ",ycurrent,"\n")
  if (control$trace>2) 
    cat("initial parameters: ",thetacurrent,"\n")
  
  ## main loop
  for (k in seq_len(control$maxit)) {
    ## get a candidate thetacand
    thetacand <- candidate.dist(thetacurrent,temps[k],control$parscale)
    ## enforce box constraints
    thetacand <- ifelse(
                        thetacand<control$lower,
                        control$lower,
                        thetacand
                        )
    thetacand <- ifelse(
                        thetacand>control$upper,
                        control$upper,
                        thetacand
                       )
    ycand <- fn(thetacand,...)/control$fnscale
    if (!is.finite(ycand)) ycand <- big
    neval <- neval+1
    
    ## see if you have a new best.params
      if (ycand<ybest) {
        ybest <- ycand
        thetabest <- thetacand
      }

    accept <- runif(1)<exp((ycurrent-ycand)/temps[k])
    if (accept) { # simulated annealing step
        thetacurrent <- thetacand
        ycurrent <- ycand
      }

    if (control$trace>1)
      cat("iter ",k," val=",ycurrent,", accept=",accept,"\n")
    if (control$trace>3) 
      cat("proposed params: ",thetacand,"\n")
    if (control$trace>2) 
      cat("current params: ",thetacurrent,"\n")

  }
    
  if (control$trace>0)
    cat("best val=",ybest,"\n")

  names(thetacurrent) <- names(thetabest) <- names(par)

  list(
       counts=c(neval,NA),
       convergence=0,
       final.params=thetacurrent,
       final.value=ycurrent,
       par=thetabest,
       value=ybest
       )
}
