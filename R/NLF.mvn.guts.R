make.lags.NLF <- function (x, lags, cov = NA, nobs = 10000) {
  x <- as.matrix(x)
  xd <- ncol(x)
  m <- length(lags)
  N <- min(nobs, nrow(x) - max(lags))
  n <- min(nobs, N)
  if (N > nobs)
    warning(" series length truncated to default in make.lags")
  start <- max(lags) + 1
  temp <- matrix(0, ncol = xd * (length(lags)), nrow = n)
  for (k in 1:length(lags)) {
    a <- start - lags[k]
    b <- a + n - 1
    temp[,(1:xd)+(k-1)*xd] <- x[(a:b),]
  }
  a <- start
  b <- a + n - 1
  if(xd == 1)
    lab <- format(paste("lag", rep(lags, rep(xd, length(lags))), sep = ""))
  else
    lab <- format(paste( rep(1:xd, length(lags)), "lag", rep(lags,rep(xd, length(lags))), sep = ""))
  dimnames(temp) <- list(NULL,lab)
  skip <- NA
  if(!is.na(cov)) {
    cov <- as.matrix(cov,drop=F)
    cov <- cov[a:b,]
    ##temp <- cbind(temp, cov[a:b,  ])
    ##cat(a, b)
    skip <- (1:ncol(cov))+m*xd
  }
  if(xd == 1)
    y <- c(x[a:b])
  else
    y <- x[a:b,]
  list(
       x=temp,
       y=y,
       nvar=m,
       cov=cov,
       lags=lags,
       skip=skip,
       start=a,
       end=b
       )
}


make.rbfbasis <- function (X, knots, fac) {
  X1 <- X-knots[1]
  nknots <- length(knots)
  if (nknots>1) {
    for (j in 2:nknots) {
      X1 <- cbind(X1,X-knots[j])
    }
  }
  exp(fac*(X1^2))
}	 


NLF.guts <- function (data.mat, model.mat, lags, nrbf = 4,
                      verbose = FALSE, plotfit = FALSE,
                      seas = FALSE, seas.data = NA, seas.sim = NA) {

  ## Version 1.0, 4 December 2007, S.P. Ellner and Bruce E. Kendall

  ## Peculiarities of the code
  ## 1. No basis functions involving cross terms
  ## 2. Model and data time series are assumed to be matrices with the same number of 
  ##    rows (= number of observation vars in data set = number of obs. vars. simulated by model)
  ## 3. Seasonal model builds an independent statistical model for each point in seasonal cycle.
  ## 4. Seasonal model not yet imported into this version 

#####################################################################################
  ## ARGUMENTS:
  ## data.mat = matrix of data time series, nobs x ntimes.data 
  ## model.mat = matrix of model time series, nobs x ntimes.sim 
  ## lags = vector of lag times for forecasting y(t) = f(y(t-lag1),y(t-lag2),....)+error
  ## nrbf = number of radial basis functions
  ## verbose: logical, print stuff out as it runs?
  ## seas: logical, is this a seasonal model?
  ## seas.data: seasonal covariate corresponding to data.matrix
  ## seas.model: seasonal covariate corresponding to data.matrix
  ##       NOTE: seas.data[j] is the covariate value used to predict data.mat[,j], so
  ##       if prediction is going to be based on a lagged covariate, this must
  ##       be set up in the calling program. 
  ## 
  ## VALUE: the NLF approximation to the negative log likelihood of the data
  ##        under the forecasting model based on model.mat. The approximation used
  ##  here is a generalized additive model for each observation variable conditional
  ##        on lagged values of all observation variables, with multivariate normal errors.     
#######################################################################################

  FAILED <- 999999; 
  if(verbose)  print("NLF multivariate version 1.0")

  multivar <- is.matrix(data.mat);  
  if (!multivar) {
    dim(data.mat) <- c(1,length(data.mat))
    dim(model.mat) <- c(1,length(model.mat))
  }
  
  nvar <- dim(data.mat)[1]
  
  ## do a lagged embedding for observation variable 1 
  data.ts <- data.mat[1,];  
  model.ts <- model.mat[1,];     

  ## Set up the RBF knots
  xm <- diff(range(model.ts))
  rbf.knots <- min(model.ts)+seq(-0.1,1.1,length=nrbf)*xm
  s <- 0.3*xm
  fac <- -1/(2*s^2)
  if (!is.finite(fac)) return(FAILED)
  if(fac==0) return(FAILED)

  ## Lag the data and set up the predicted values & seasonal indices
  Lags.model <- make.lags.NLF(model.ts,lags=lags,cov=seas.sim)
  Lags.data <- make.lags.NLF(data.ts,lags=lags,cov=seas.data)

  data.pred <- matrix(Lags.data$y,ncol=1)
  model.pred <- matrix(Lags.model$y,ncol=1)

  if (verbose) {
    print("calculating ridge functions")
    print(date())
  }

  rbfbasis.model <- make.rbfbasis(Lags.model$x,knots=rbf.knots,fac=fac)
  rbfbasis.data <- make.rbfbasis(Lags.data$x,knots=rbf.knots,fac=fac)

  if (multivar) {
    for (jvar in 2:nvar) {
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

  prediction.errors <- matrix(0,dim(data.pred)[1],nvar)
  model.residuals <- matrix(0,dim(model.pred)[1],nvar)
  for (jvar in 1:nvar)  {
    model.lm <- lm(model.pred[,jvar]~rbfbasis.model-1)
    model.residuals[,jvar] <- residuals(model.lm)
    ck <- as.vector(coef(model.lm))
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
    LQL <- sum(dnorm(prediction.errors[,1],mean=0,sd=sigma.model,log=TRUE),na.rm=T)
  } else {
    sigma.model <- cov(model.residuals)
    ## the following does not require package 'mvtnorm'
    distval <- mahalanobis(
                           x=prediction.errors,
                           center=0,
                           cov=sigma.model
                           )
    logdet <- sum(log(eigen(sigma.model,symmetric=TRUE,only.values=TRUE)$values))
    LQL <- -0.5*sum(nvar*log(2*pi)+logdet+distval,na.rm=TRUE)
    ##  the following requires the package 'mvtnorm'
    ##    LQL <- sum(dmvnorm(prediction.errors,sigma=sigma.model,log=TRUE),na.rm=T)
  }  

  -LQL
}
