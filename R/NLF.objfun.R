NLF.objfun <- function (params.fitted, object, params, par.index, ts.index, xstart, 
                        times, lags, seed, nrbf = 4, verbose = FALSE) {
  
  FAILED <- 999999
  params[par.index] <- params.fitted
  params <- as.matrix(params)
  xstart <- as.matrix(xstart)
  
  ## Need to extract number of state variables (nvar) from pomp object
  ## Need to include simulation times in problem specification

  ## Evaluates the NLF objective function given a POMP object.

  ## Version 0.1, 3 Dec. 2007, Bruce E. Kendall & Steve P. Ellner

  data.ts <- as.numeric(data.array(object)[ts.index,])
  
  ## set the random seed (be very careful about this)
  if (!is.null(seed)) {
    if (!exists('.Random.seed',envir=.GlobalEnv)) runif(1)
    save.seed <- get('.Random.seed',envir=.GlobalEnv)
    set.seed(seed)
  }

  x <- rprocess(object,xstart=xstart,times=times,params=params) # simulate the process model
  if (!all(is.finite(x))) return(FAILED)

  y <- rmeasure(object,x=x[,,-1,drop=F],times=times[-1],params=params)
  
  if (!is.null(seed)) assign('.Random.seed',save.seed,envir=.GlobalEnv)  # restore the RNG state
  
  ## Test whether the model time series is valid
  
  if (all(is.finite(y))) {
    model.ts <- as.numeric(y[ts.index,1,])
#    if (verbose) print(summary(model.ts))
    
    LQL <- NLF.guts(data.mat=data.ts, model.mat=model.ts, lags=lags, nrbf=nrbf, verbose=F, plotfit=F)
  } else {
    LQL <- FAILED
  }
  if (verbose) print(c(LQL, params))

  LQL
}
