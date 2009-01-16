NLF.LQL <- function (params.fitted, object, params, par.index,
                     times, lags, period, tensor, seed, nrbf = 4, verbose = FALSE,
                     bootstrap = FALSE, bootsamp = NULL) {

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Computes the vector of log quasi-likelihood values at the observations  
# Note that the log QL itself is returned, not the negative log QL,
# so a large NEGATIVE value is used to flag bad parameters 
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  FAILED =  -99999999999
  params[par.index] <- params.fitted
  
  params <- as.matrix(params)
  xstart <- init.state(object,params)

  ## Need to extract number of state variables (nvar) from pomp object
  ## Need to include simulation times in problem specification
  ## Evaluates the NLF objective function given a POMP object.
  ## Version 0.1, 3 Dec. 2007, Bruce E. Kendall & Stephen P. Ellner
  ## Version 0.2, May 2008, Stephen P. Ellner  

  data.ts <- data.array(object)

  ## set the random seed (be very careful about this)
  if (!is.null(seed)) {
    if (!exists('.Random.seed',envir=.GlobalEnv)) runif(1)
    save.seed <- get('.Random.seed',envir=.GlobalEnv)
   set.seed(seed)
  }

  x <- try(
           rprocess(object,xstart=xstart,times=times,params=params), # simulate the process model
           silent=FALSE
           )
  if (inherits(x,'try-error'))
    stop("'NLF.LQL' reports: error in user 'rprocess'")

  if (!all(is.finite(x))) return(FAILED)

  y <- try(
           rmeasure(object,x=x[,,-1,drop=F],times=times[-1],params=params),
           silent=FALSE
           )
  if (inherits(y,'try-error'))
    stop("'NLF.LQL' reports: error in user 'rmeasure'")

  if (!is.null(seed)) assign('.Random.seed',save.seed,envir=.GlobalEnv)  # restore the RNG state

 ## Test whether the model time series is valid
  model.ts <- array(dim=c(nrow(data.ts),length(times)-1),dimnames=list(rownames(data.ts),NULL))
  if (!all(is.finite(x))) return(FAILED)

  model.ts[,] <- y[,1,]
  LQL <- try(
             NLF.guts(
                      data.mat=data.ts,
                      data.times=time(object),
                      model.mat=model.ts,
                      model.times=times[-1],
                      lags=lags,
                      period=period, tensor=tensor, 
                      nrbf=nrbf,
                      verbose=F,
                      bootstrap,
                      bootsamp,
                      plotfit=F
                      ),
             silent=FALSE
             )
  if (inherits(LQL,"try-error"))
    stop("'NLF.LQL' reports: error in 'NLF.guts'")
  LQL 
 }

