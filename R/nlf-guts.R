NLF.guts <- function (data.mat, data.times, model.mat, model.times, lags, period,
                      tensor, nrbf = 4, verbose = FALSE, plotfit = FALSE,
                      bootstrap = FALSE, bootsamp = NULL) {

  ## Version 1.0, 4 December 2007, S.P. Ellner and Bruce E. Kendall
  ## Verstion 1.1, 19 June 2008, S.P. Ellner and A.A. King. 

  ## Peculiarities of the code
  ## 1. No basis functions involving cross terms of state variables 
  ## 2. Model and data time series are assumed to be matrices with the same number of 
  ##    rows (= number of observation vars in data set = number of obs. vars. simulated by model)
  ##    This is formally a requirement of pomp objects, and here it is absolutely required.
  

#####################################################################################
  ## ARGUMENTS:
  ## data.mat = matrix of data time series, nobs x ntimes.data 
  ## model.mat = matrix of model time series, nobs x ntimes.sim 
  ## lags = vector of lag times for forecasting y(t) = f(y(t-lag1),y(t-lag2),....)+error
  ## nrbf = number of radial basis functions
  ## verbose: logical, print stuff out as it runs?
  ## period: period=NA means the model is nonseasonal. period=integer>0 is the period of
  ##         seasonal forcing in 'real time' (the units of model.times). 
  ## tensor: logical. if FALSE, the fitted model is a gam with time(mod period) as
  ##         one of the predictors, i.e. a gam with time-varying intercept. 
  ##         if TRUE, the fitted model is a gam with lagged state variables as
  ##         predictors and time-periodic coefficients, constructed using tensor
  ##         products of basis functions of state variables with basis functions of time. 
  ##
  ##       NOTE: periodic models are constructed so that the time variable on the
  ##       right hand side corresponds to the observation time of the predictee,
  ##       not of the predictors 
  ##
  ## 
  ## VALUE: the NLF approximation to the log likelihood of the data series
  ##        under the forecasting model based on model.mat. The approximation used
  ##        here is a generalized additive model for each observation variable, conditional
  ##        on lagged values of all observation variables, with multivariate normal errors.     
  ##        The return from this function is the vector of log quasi-likelihood values at 
  ##        the data points; this must be summed to get the log quasiliklihood function.  
  ##
  ## IMPORTANT NOTE TO FUTURE PROGRAMMERS: It may appear at first glance that basis
  ## functions for the data series, and other things related to the data series, could be
  ## constructed once and for all rather than rebuilt on each call to this function. 
  ## THIS IS NOT TRUE. The basis functions are constructed anew on each call using
  ## information from the model-simulated time series, and this feature is important 
  ## for reliable NLF parameter estimation because it rules out spurious good fits
  ## that really are due to extrapolation far from the model-simulated time series.   
#######################################################################################


  FAILED = -999999999; 
  if (verbose) print("NLF multivariate version 1.0")

  nvar <- nrow(data.mat)
  multivar <- (nvar>1)

  ## do a lagged embedding for observation variable 1 
  data.ts <- data.mat[1,]
  model.ts <- model.mat[1,]

  ## Set up the RBF knots
  xm <- diff(range(model.ts))
  rbf.knots <- min(model.ts)+seq(-0.1,1.1,length=nrbf)*xm
  s <- 0.3*xm
  fac <- -1/(2*s^2)
  if (!is.finite(fac)) return(FAILED)
  if (fac==0) return(FAILED)

  seas <- (!is.na(period)&&abs(period)>0)

  if (seas) {
    seas.sim <- model.times%%abs(period)
    seas.data <- data.times%%abs(period)
  } else {
    seas.sim <- NULL
    seas.data <- NULL
  }

  ## Lag the data and set up the predicted values & seasonal indices
  Lags.model <- make.lags.NLF(model.ts,lags=lags,cov=seas.sim)
  Lags.data <- make.lags.NLF(data.ts,lags=lags,cov=seas.data)

  if (bootstrap) {
    Lags.data$x <- Lags.data$x[bootsamp,]
    Lags.data$y <- Lags.data$y[bootsamp]
    if (seas) 
      Lags.data$cov <- Lags.data$cov[bootsamp]
  }	

  data.pred <- matrix(Lags.data$y,ncol=1)
  model.pred <- matrix(Lags.model$y,ncol=1)

  if (verbose) {
    print("calculating ridge functions")
    print(date())
  }

  rbfbasis.model <- make.rbfbasis(Lags.model$x,knots=rbf.knots,fac=fac)
  rbfbasis.data <- make.rbfbasis(Lags.data$x,knots=rbf.knots,fac=fac)

  if (seas) {
    ## Set up the RBF knots
    rbf.cov.knots <- seq(-0.1,1.1,length=nrbf)*period
    s <- 0.3*period
    fac.cov <- -1/(2*s^2)
    if (!is.finite(fac.cov)) return(FAILED)
    if (fac.cov==0) return(FAILED)
    
    rbfbasis.cov.model <- make.rbfbasis(Lags.model$cov,knots=rbf.cov.knots,fac=fac.cov)
    rbfbasis.cov.data <- make.rbfbasis(Lags.data$cov,knots=rbf.cov.knots,fac=fac.cov)
  }

  if (multivar) {
    for (jvar in seq(from=2,to=nvar,by=1)) {
      data.ts <- data.mat[jvar,]
      model.ts <- model.mat[jvar,]

      ## Set up the RBF knots
      xm <- diff(range(model.ts))
      rbf.knots <- min(model.ts)+seq(-0.1,1.1,length=nrbf)*xm
      s <- 0.3*xm
      fac <- -1/(2*s^2)
      if (fac==0) return(FAILED)
   
      ## Lag the data and set up the predicted values & seasonal indices
      Lags.model <- make.lags.NLF(model.ts,lags=lags,cov=seas.sim)
      Lags.data <- make.lags.NLF(data.ts,lags=lags,cov=seas.data)
  
      if (bootstrap) {
	Lags.data$x=Lags.data$x[bootsamp,]
	Lags.data$y=Lags.data$y[bootsamp]
      }	

      data.pred <- cbind(data.pred,Lags.data$y)
      model.pred <- cbind(model.pred,Lags.model$y)
      rbfbasis.model <- cbind(rbfbasis.model,make.rbfbasis(Lags.model$x,knots=rbf.knots,fac=fac))
      rbfbasis.data <- cbind(rbfbasis.data,make.rbfbasis(Lags.data$x,knots=rbf.knots,fac=fac))
 
      if (verbose) {
        print("done with ridge functions")
        print(date())
      }
    }
  }

  if (seas) {
      if(tensor) { 
         # make gam coefficients time-dependent
         rbfbasis.model <- make.tensorbasis.NLF(rbfbasis.model,rbfbasis.cov.model)      
	 rbfbasis.data <- make.tensorbasis.NLF(rbfbasis.data,rbfbasis.cov.data)
      }else{ 
         # add time-varying intercept  
         rbfbasis.model <- cbind(rbfbasis.model,rbfbasis.cov.model)
         rbfbasis.data <- cbind(rbfbasis.data,rbfbasis.cov.data)
      }
  }
  
  prediction.errors <- matrix(0,dim(data.pred)[1],nvar)
  model.residuals <- matrix(0,dim(model.pred)[1],nvar)

  for (jvar in seq_len(nvar)) {
    model.lm <- .lm.fit(rbfbasis.model,model.pred[,jvar])
    model.residuals[,jvar] <- model.lm$residuals
    ck <- model.lm$coefficients
    if (verbose) {
      print(ck)
      print(summary(model.lm))
    }  

    fitted.data <- rbfbasis.data%*%matrix(ck,ncol=1)
    prediction.errors[,jvar] <- data.pred[,jvar]-fitted.data
    if (plotfit) {
      par(mfrow=c(1,2))
      plot(model.pred,fitted(model.lm),pch='.')
      abline(0,1)
      plot(data.pred,fitted.data)
      abline(0,1)
      par(mfrow=c(1,1))
    }
  }

  if (nvar==1) {
    sigma.model <- sd(model.residuals[,1])
    LQL <- dnorm(prediction.errors[,1],mean=0,sd=sigma.model,log=TRUE)
  } else {
    sigma.model <- cov(model.residuals)
    LQL <- dmvnorm(prediction.errors,sigma=sigma.model,log=TRUE) ## NOTE: This could be improved using GLS.
  }  
  
  LQL
}
