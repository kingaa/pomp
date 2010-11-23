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

  npar <- length(par)
  neval <- 0
  
  control.default <- list(
                          maxit=10000,
                          temp=1,
                          tmax=10,
                          sched=NULL,
                          candidate.dist=NULL,
                          parscale=1,
                          lower=-Inf,
                          upper=Inf,
                          trace=0
                          )
  control.default[names(control)] <- control
  control <- control.default

  if (is.null(control$sched))           # default cooling schedule
    control$sched <- function (k) control$temp/log(((k-1)%/%control$tmax)*control$tmax+exp(1))

  if (is.function(control$sched))
    temps <- sapply(seq.int(from=1,to=control$maxit,by=1),control$sched)
  else if (is.numeric(control$sched))
    temps <- control$sched
  
  if (is.null(control$candidate.dist))
    candidate.dist <- function (temp) rnorm(n=npar,mean=0,sd=control$parscale*temp)

  if (length(control$lower)<npar)
    control$lower <- rep(control$lower,npar)
  if (length(control$upper)<npar)
    control$upper <- rep(control$upper,npar)

  ## initialization for the algorithm
  laststep <- 0
  thetabest <- thetacurrent <- par
  ybest <- ycurrent <- fn(thetacurrent,...)
  neval <- 1
  
  if (control$trace>0)
    cat("initial evaluation: ",ycurrent,"\n")
  if (control$trace>2) 
    cat("initial parameters: ",thetacurrent,"\n")
  
  ## main loop
  for (k in seq_len(control$maxit)) {
    ## get a candidate thetacand
    thetacand <- thetacurrent+candidate.dist(temps[k])
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
    ycand <- fn(thetacand,...)
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

  list(
       counts=c(neval,NA),
       convergence=0,
       final.params=thetacurrent,
       final.value=ycurrent,
       par=thetabest,
       value=ybest
       )
}
