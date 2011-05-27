nlf <- function (object, start, est, lags,
                 period = NA, tensor = FALSE,
                 nconverge = 1000, nasymp = 1000, 
                 seed = 1066, transform = function(x)x,
                 nrbf = 4, method = "subplex",
                 skip.se = FALSE, verbose = FALSE, gr = NULL, 
                 bootstrap = FALSE, bootsamp = NULL,
                 lql.frac = 0.1, se.par.frac = 0.1,
                 eval.only = FALSE, ...) {

  ## Fit a POMP object using NLF
  ## v. 0.1, 3 Dec. 2007 
  ## by Bruce Kendall & Steve Ellner
  ## 
  ## v. 0.2, 30 May 2008, by Steve Ellner 
  ## Adds automatic Wald asymptotic standard errors and the 
  ## capability for moving-blocks bootstrap standard errors.  
  ## Quadratic regression near optimum used to select increments
  ## for finite-difference approximations to gradient and Hessian 
  ##
  ## v 1.0, 19 June 2008 by Steve Ellner and Aaron King
  ## adds capacity to fit models with periodically time-varying parameters
  ## of known period and improves the compatibility with the standard for pomp objects

  if (!is(object,'pomp'))
    stop("'object' must be a 'pomp' object")

  if (!is.function(transform))
    stop(sQuote("transform")," must be a function!")

  if (eval.only) est <- 1

  if (is.character(est)) {
    if (!all(est%in%names(start)))
      stop(sQuote("nlf")," error: parameters named in ",sQuote("est")," must exist in ",sQuote("start"))

    par.index <- which(names(start)%in%est)

  } else if (is.numeric(est)) {
    est <- as.integer(est)
    if (any((est<1)|(est>length(start))))
      stop(sQuote("nlf")," error: trying to estimate parameters that don't exist!")
    par.index <- as.integer(est)
  }

  if (!(is.numeric(lql.frac)&&(lql.frac>0)&&(lql.frac<1)))
    stop(sQuote("nlf")," error: ",sQuote("lql.frac")," must be in (0,1)")

  if (!(is.numeric(se.par.frac)&&(se.par.frac>0)&&(se.par.frac<1)))
    stop(sQuote("nlf")," error: ",sQuote("se.par.frac")," must be in (0,1)")

  params <- start
  guess <- params[par.index]

  dt.tol <- 1e-3
  times <- time(object,t0=FALSE)
  t0 <- timezero(object)
  dt <- diff(times)
  if (diff(range(dt))>dt.tol*mean(dt))
    stop(sQuote("nlf")," requires evenly spaced sampling times")
  dt <- times[2]-times[1]

  ## Vector of times to output the simulation
  times <- seq(
               from=t0+nconverge*dt,
               length=nasymp,
               by=dt
               )

  if (eval.only) {
    val <- nlf.objfun(
                      params.fitted=guess,
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
                      )
    return(-val)
  }

  if (method == 'subplex') {
    opt <- subplex::subplex(
                            par=guess,
                            fn=nlf.objfun,
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
                            bootsamp=bootsamp,
                            control=list(...)
                            )
  } else {
    opt <- optim(
                 par=guess,
                 fn=nlf.objfun,
                 gr=gr,
                 method=method, 
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
                 bootsamp=bootsamp,
                 control=list(...)
                 )  
  }

  opt$est <- est
  opt$value <- -opt$value
  params[par.index] <- opt$par
  opt$params <- params
  opt$par <- NULL

  if (!skip.se) { ## compute estimated Variance-Covariance matrix of fitted parameters 
    fitted <- params[par.index]
    nfitted <- length(fitted)
    Jhat <- matrix(0,nfitted,nfitted)
    Ihat <- Jhat
    f0 <- NLF.LQL(fitted,object=object, params=params, par.index=par.index, 
                  times=times, t0=t0,
                  lags=lags, period=period, tensor=tensor, seed=seed,
                  transform=transform, nrbf=4, 
                  verbose=FALSE)
    F0 <- mean(f0,na.rm=T)

    npts <- length(f0)
    nlags <- round(5*npts^0.2) ## Number of lags to use in Newey-West covariance estimator 

    ## find a good epsilon 
    h <- se.par.frac
    if (verbose)
      cat("h in NLF = ", h, "\n")
    eps <- rep(h,nfitted)

    for (i in seq_len(nfitted)) {
      Fvals <- rep(0,5)
      Fvals[3] <- F0  
      guess <- fitted
      guess[i] <- fitted[i]-sqrt(2)*h*abs(fitted[i])  
      Fvals[1] <- mean(NLF.LQL(guess,object=object, params=params, par.index=par.index, 
                               times=times, t0=t0, lags=lags, period=period, tensor=tensor,
                               seed=seed, transform=transform,
                               nrbf=4, verbose=FALSE),na.rm=T)
      guess <- fitted
      guess[i] <- fitted[i]-h*abs(fitted[i])
      Fvals[2] <- mean(NLF.LQL(guess,object=object, params=params, par.index=par.index, 
                               times=times, t0=t0, lags=lags, period=period, tensor=tensor,
                               seed=seed, transform=transform, nrbf=4, 
                               verbose=FALSE),na.rm=T)
      guess <- fitted
      guess[i] <- fitted[i]+h*abs(fitted[i])
      Fvals[4] <- mean(NLF.LQL(guess,object=object, params=params, par.index=par.index, 
                               times=times, t0=t0, lags=lags, period=period, tensor=tensor,
                               seed=seed, transform=transform, nrbf=4, 
                               verbose=FALSE),na.rm=T)
      guess <- fitted
      guess[i] <- fitted[i]+sqrt(2)*h*abs(fitted[i])
      Fvals[5] <- mean(NLF.LQL(guess,object=object, params=params, par.index=par.index, 
                               times=times, t0=t0, lags=lags, period=period, tensor=tensor,
                               seed=seed, transform=transform, nrbf=4, 
                               verbose=FALSE),na.rm=T)
      FAILED =  - 999999
      Fvals[Fvals < FAILED+10] <- NA
      xvals <- c(sqrt(2),1,0,1,sqrt(2))*h*fitted[i]
      c2 <- lm(Fvals~I(xvals^2))$coef[2]
      eps[i] <- sqrt(abs(lql.frac/c2))
    }

    if (verbose)
      cat("epsilon in NLF =",t(eps), "\n")

    Imat <- matrix(0,npts,nfitted)
    for (i in seq_len(nfitted)) {
      guess.up <- fitted
      guess.up[i] <- guess.up[i]+eps[i]
      f.up <- NLF.LQL(guess.up,object=object, params=params, par.index=par.index, 
                      times=times, t0=t0, lags=lags, period=period, tensor=tensor,
                      seed=seed, transform=transform, nrbf=4, 
                      verbose=FALSE)
      F.up <- mean(f.up,na.rm=T)

      f.up2 <- NLF.LQL(guess.up,object=object, params=params, par.index=par.index, 
                       times=times, t0=t0, lags=lags, period=period, tensor=tensor,
                       seed=seed, transform=transform, nrbf=4, 
                       verbose=FALSE)

      if (verbose)
        cat("Fitted param ", i, F.up, mean(f.up2,na.rm=T)," up in ",sQuote("nlf"),"\n")

      guess.down <- fitted
      guess.down[i] <- guess.down[i]-eps[i]
      f.down <- NLF.LQL(guess.down,object=object, params=params, par.index=par.index, 
                        times=times, t0=t0, lags=lags, period=period, tensor=tensor,
                        seed=seed, transform=transform, nrbf=4, 
                        verbose=FALSE)
      F.down <- mean(f.down,na.rm=T)

      if (verbose)
        cat("Fitted param ",i, F.down," down in ",sQuote("NLF"),"\n")

      Jhat[i,i] <- (F.up + F.down-2*F0)/(eps[i]*eps[i])
      Imat[,i] <- (f.up-f.down)/(2*eps[i])
      Ihat[i,i] <- Newey.West(Imat[,i],Imat[,i],nlags)
    }

    for (i in seq_len(nfitted-1)) {
      for (j in seq(from=i+1,to=nfitted,by=1)) {
        guess.uu <- fitted
        guess.uu[i] <- guess.uu[i]+eps[i]
        guess.uu[j] <- guess.uu[j]+eps[j]
        F.uu <- mean(NLF.LQL(guess.uu,object=object, params=params, par.index=par.index,
                             times=times, t0=t0, lags=lags, period=period, tensor=tensor,
                             seed=seed, transform=transform, nrbf=4, 
                             verbose=FALSE),na.rm=T)

        guess.ud <- fitted 
        guess.ud[i] <- guess.ud[i]+eps[i]
        guess.ud[j] <- guess.ud[j]-eps[j]
        F.ud <- mean(NLF.LQL(guess.ud,object=object, params=params, par.index=par.index,
                             times=times, t0=t0, lags=lags, period=period, tensor=tensor,
                             seed=seed, transform=transform, nrbf=4, 
                             verbose=FALSE),na.rm=T) 

        guess.du <- fitted 
        guess.du[i] <- guess.du[i]-eps[i]
        guess.du[j] <- guess.du[j]+eps[j]
        F.du <- mean(NLF.LQL(guess.du,object=object, params=params, par.index=par.index,
                             times=times, t0=t0, lags=lags, period=period, tensor=tensor,
                             seed=seed, transform=transform, nrbf=4, 
                             verbose=FALSE),na.rm=T) 

        guess.dd <- fitted 
        guess.dd[i] <- guess.dd[i]-eps[i]
        guess.dd[j] <- guess.dd[j]-eps[j] 
        F.dd <- mean(NLF.LQL(guess.dd,object=object, params=params, par.index=par.index,
                             times=times, t0=t0, lags=lags, period=period, tensor=tensor,
                             seed=seed, transform=transform, nrbf=4,
                             verbose=FALSE),na.rm=T) 

        dij <- (F.uu+F.dd)-(F.ud+F.du)
        dij <- dij/(4*eps[i]*eps[j]) 
        Jhat[i,j] <- dij
        Jhat[j,i] <- dij
        Ihat[i,j] <- Newey.West(Imat[,i],Imat[,j],nlags)
        Ihat[j,i] <- Ihat[i,j]  
      }
    }
    opt$Jhat <- Jhat
    opt$Ihat <- Ihat
    negJinv <- -solve(Jhat)
    Qhat <- negJinv%*%Ihat%*%negJinv
    opt$Qhat <- Qhat
    opt$se <- sqrt(diag(Qhat))/sqrt(npts)
    names(opt$se) <- names(start)[par.index]
    opt$npts <- npts
  }
  
  opt
}

