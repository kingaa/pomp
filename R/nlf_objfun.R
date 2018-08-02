## Evaluates the NLF objective function given a POMP object.
## Version 0.1, 3 Dec. 2007, Bruce E. Kendall & Stephen P. Ellner
## Version 0.2, May 2008, Stephen P. Ellner  
###
### Computes the vector of log quasi-likelihood values at the observations  
### Note that the log QL itself is returned, not the negative log QL,
### so a large NEGATIVE value is used to flag bad parameters 
###

nlf.objfun <- function (...) 
    -sum(nlf.lql(...),na.rm=TRUE)

nlf.lql <- function (params.fitted, object, params, par.index,
                     transform = FALSE,
                     times, t0, lags, period, tensor, seed = NULL,
                     transform.data = identity, nrbf = 4, verbose = FALSE,
                     bootstrap = FALSE, bootsamp = NULL) {

    ep <- paste0("in ",sQuote("nlf.lql"),": ")
    
    transform <- as.logical(transform)

    FAILED =  -99999999999
    params[par.index] <- params.fitted
    
    if (transform)
        params <- partrans(object,params,dir="fromEst")

    data.ts <- obs(object)
    
    y <- tryCatch(
        simulate(object,times=times,t0=t0,params=params,seed=seed,
                 obs=TRUE,states=FALSE),
        error = function (e) {
            stop(ep,"simulation error: ",conditionMessage(e),call.=FALSE) # nocov
        }
    )

    ## Test whether the model time series is valid
    if (!all(is.finite(y))) return(FAILED)

    model.ts <- array(
        dim=c(nrow(data.ts),length(times)),
        dimnames=list(rownames(data.ts),NULL)
    )
    model.ts[,] <- apply(y[,1,,drop=FALSE],c(2,3),transform.data)
    data.ts[,] <- apply(data.ts,2,transform.data)
    
    tryCatch(
        nlf.guts(
            data.mat=data.ts,
            data.times=time(object),
            model.mat=model.ts,
            model.times=times,
            lags=lags,
            period=period,
            tensor=tensor, 
            nrbf=nrbf,
            bootstrap=bootstrap,
            bootsamp=bootsamp
        ),
        error = function (e) {
            stop(ep,conditionMessage(e),call.=FALSE) # nocov
        }
    )
}
