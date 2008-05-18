NLF <- function(object, xstart, params, est, xfit, lags, nconverge = 1000, nasymp = 1000, 
                seed = 1066, nrbf = 4, method = 'subplex', verbose = FALSE, gr = NULL, ...) {

  ## Fit a POMP object using NLF
  ## v. 0.1, 3 Dec. 2007 
  ## by Bruce Kendall & Steve Ellner

## >>>>>>>>> Stolen from trajmatch.R <<<<<<<<<<<<<<<<<<<<<<
  if (!inherits(object,'pomp'))
    stop("'object' must be a 'pomp' object")
  if (!is.character(est))
    stop("'est' must be a character vector containing the names of the parameters to be estimated")
  par.est <- which(names(params)%in%est)
  par.index <- seq(along=par.est)
  guess <- params[par.est]
# >>>>>>>>> End Stolen from trajmatch.R <<<<<<<<<<<<<<<<<<<<<<
  if (!is.character(xfit))
    stop("'xfit' must be a character vector containing the names of the state variables to be fit")
  x.xfit <- which(names(xstart)%in%xfit)
  ts.index <- seq(along=x.xfit)
  
  
  times <- c(0,nconverge+(1:nasymp)) # Vector of times to output the simulation
  
  
  if (method == 'subplex') {
    opt <- subplex(guess, NLF.objfun, object=object, params=params, par.index=par.index, 
                   ts.index=ts.index, xstart=xstart, times=times, lags=lags, seed=seed, nrbf=nrbf, 
                   verbose=verbose)
  } else {
    opt <- optim(par=guess, fn=NLF.objfun, gr=gr, method=method, control=list(...), 
                 object=object, params=params, par.index=par.index, 
                 ts.index=ts.index, xstart=xstart, times=times, lags=lags, seed=seed, nrbf=nrbf, 
                 verbose=verbose)  
  }

#print(NLF.objfun(guess,object=object, params=params, par.index=par.index, 
#      ts.index=ts.index, xstart=xstart, times=times, lags=lags, seed=seed, nrbf=4, 
#      verbose=verbose))

# >>>>>>>>> Stolen from trajmatch.R <<<<<<<<<<<<<<<<<<<<<<  
  opt$value <- -opt$value
  params[par.est] <- opt$par[par.index]
  opt$xstart <- xstart
  opt$params <- params
  opt$par <- NULL
  opt
# >>>>>>>>> End Stolen from trajmatch.R <<<<<<<<<<<<<<<<<<<<<<
}
