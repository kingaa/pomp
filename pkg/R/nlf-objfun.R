nlf.objfun <- function (params.fitted, object, params, par.index,
                        times, lags, period, tensor, seed, transform = function(x)x,
                        nrbf = 4, verbose = FALSE, bootstrap = FALSE, bootsamp = NULL) 
{
  -sum(
       NLF.LQL(
               params.fitted=params.fitted,
               object=object,
               params=params,
               par.index=par.index,
               times=times,
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
