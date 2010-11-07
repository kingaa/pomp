NLF.LQL <- function (params.fitted, object, params, par.index,
                     times, t0, lags, period, tensor, seed = NULL, transform = identity,
                     nrbf = 4, verbose = FALSE,
                     bootstrap = FALSE, bootsamp = NULL) {

###>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
### Computes the vector of log quasi-likelihood values at the observations  
### Note that the log QL itself is returned, not the negative log QL,
### so a large NEGATIVE value is used to flag bad parameters 
###>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  FAILED =  -99999999999
  params[par.index] <- params.fitted
  
  params <- as.matrix(params)

  ## Need to extract number of state variables (nvar) from pomp object
  ## Need to include simulation times in problem specification
  ## Evaluates the NLF objective function given a POMP object.
  ## Version 0.1, 3 Dec. 2007, Bruce E. Kendall & Stephen P. Ellner
  ## Version 0.2, May 2008, Stephen P. Ellner  

  data.ts <- obs(object)
  
  y <- try(
           simulate(object,times=times,t0=t0,params=params,seed=seed,obs=TRUE,states=FALSE),
           silent=FALSE
           )
  if (inherits(y,"try-error"))
    stop(sQuote("NLF.LQL")," reports: error in simulation")
  ## Test whether the model time series is valid
  if (!all(is.finite(y))) return(FAILED)

  model.ts <- array(
                    dim=c(nrow(data.ts),length(times)),
                    dimnames=list(rownames(data.ts),NULL)
                    )
  model.ts[,] <- apply(y[,1,,drop=FALSE],c(2,3),transform)
  data.ts[,] <- apply(data.ts,2,transform)
  
  LQL <- try(
             NLF.guts(
                      data.mat=data.ts,
                      data.times=time(object),
                      model.mat=model.ts,
                      model.times=times,
                      lags=lags,
                      period=period,
                      tensor=tensor, 
                      nrbf=nrbf,
                      verbose=FALSE,
                      bootstrap,
                      bootsamp,
                      plotfit=FALSE
                      ),
             silent=FALSE
             )
  if (inherits(LQL,"try-error"))
    stop(sQuote("NLF.LQL")," reports: error in ",sQuote("NLF.guts"))
  LQL 
}



nlf.objfun <- function (params.fitted, object, params, par.index,
                        times, t0, lags, period, tensor, seed, transform = function(x)x,
                        nrbf = 4, verbose = FALSE, bootstrap = FALSE, bootsamp = NULL) 
{
  -sum(
       NLF.LQL(
               params.fitted=params.fitted,
               object=object,
               params=params,
               par.index=par.index,
               times=times,
               t0=t0,
               lags=lags,
               period=period,
               tensor=tensor,
               seed=seed,
               transform=transform,
               nrbf=nrbf,
               verbose=verbose,
               bootstrap=bootstrap,
               bootsamp=bootsamp
               ),
       na.rm=T
       )
} 
