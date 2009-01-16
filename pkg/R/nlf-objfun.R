nlf.objfun <- function (params.fitted, object, params, par.index,
                        times, lags, period, tensor, seed, nrbf = 4,
                        verbose = FALSE, bootstrap = FALSE, bootsamp = NULL) 
{
  -sum(
       NLF.LQL(
               params.fitted,
               object,
               params,
               par.index,
               times,
               lags,
               period, tensor,
               seed,
               nrbf,
               verbose,
               bootstrap,
               bootsamp
               ),
       na.rm=T
       )
} 
