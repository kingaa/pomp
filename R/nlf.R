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

setClass(
  "nlfd_pomp",
  contains="pomp",
  slots=c(
    transform = "logical",
    transform.data = "function",
    est = 'character',
    lags="integer",
    nconverge = 'integer',
    nasymp = 'integer',
    seed="integer",
    period="numeric",
    tensor="logical",
    nrbf="integer",
    method="character",
    lql.frac="numeric",
    se.par.frac="numeric",
    Qhat="matrix",
    se="numeric",
    logql="numeric"
  ),
  prototype=prototype(
    transform=FALSE,
    transform.data=identity,
    est=character(0),
    lags=integer(0),
    nconverge=0L,
    nasymp=0L,
    seed=0L,
    period=as.numeric(NA),
    tensor=FALSE,
    nrbf=4L,
    method=character(0),
    lql.frac=0.1,
    se.par.frac=0.1,
    Qhat=matrix(NA,0,0),
    se=numeric(0),
    logql=as.numeric(NA)
  )
)

setMethod(
  "nlf",
  signature=signature(object="pomp"),
  definition=function (object,
    start, est, lags,
    period = NA, tensor = FALSE,
    nconverge = 1000L, nasymp = 1000L,
    seed = 1066L, transform.data,
    nrbf = 4L,
    method = c(
      "subplex", "Nelder-Mead", "BFGS", "CG",
      "L-BFGS-B", "SANN", "Brent"
    ),
    skip.se = FALSE,
    verbose = getOption("verbose"),
    bootsamp = NULL,
    lql.frac = 0.1, se.par.frac = 0.1,
    eval.only = FALSE,
    transform = FALSE, ...)
  {
    ep <- paste0("in ",sQuote("nlf"),": ")
    transform <- as.logical(transform)
    if (missing(transform.data)) transform.data <- identity
    transform.data <- match.fun(transform.data)
    period <- as.numeric(period)
    tensor <- as.logical(tensor)
    skip.se <- as.logical(skip.se)
    eval.only <- as.logical(eval.only)
    seed <- as.integer(seed)
    lql.frac <- as.numeric(lql.frac)
    se.par.frac <- as.numeric(se.par.frac)
    bootstrap <- !is.null(bootsamp)
    bootsamp <- as.integer(bootsamp)
    lags <- as.integer(lags)
    nrbf <- as.integer(nrbf)
    nasymp <- as.integer(nasymp)
    nconverge <- as.integer(nconverge)

    method <- match.arg(method)

    if (missing(est) || length(est)==0) eval.only <- TRUE
    if (eval.only) est <- character(0)
    if (missing(start)) start <- coef(object)
    if (is.list(start)) start <- unlist(start)
    if (!is.character(est))
      stop(ep,sQuote("est")," must name the parameters to be estimated",
        call.=FALSE)
    if (!all(est%in%names(start)))
      stop(ep,"parameters named in ",sQuote("est"),
        " must exist in ",sQuote("start"),call.=FALSE)

    nlf.internal(
      object=object,
      start=start,
      est=est,
      lags=lags,
      period=period,
      tensor=tensor,
      nconverge=nconverge,
      nasymp=nasymp,
      seed=seed,
      nrbf=nrbf,
      method=method,
      skip.se=skip.se,
      verbose=verbose,
      bootstrap=bootstrap,
      bootsamp=bootsamp,
      lql.frac=lql.frac,
      se.par.frac=se.par.frac,
      eval.only=eval.only,
      transform=transform,
      transform.data=transform.data,
      ...
    )
  }
)

setMethod(
  "nlf",
  signature=signature(object="nlfd_pomp"),
  definition=function (object, start, est, lags,
    period, tensor, nconverge, nasymp, seed,
    transform.data, nrbf, method, lql.frac, se.par.frac,
    transform, ...)
  {
    if (missing(start)) start <- coef(object)
    if (is.list(start)) start <- unlist(start)
    if (missing(est)) est <- object@est
    if (missing(lags)) lags <- object@lags
    if (missing(period)) period <- object@period
    if (missing(tensor)) tensor <- object@tensor
    if (missing(nconverge)) nconverge <- object@nconverge
    if (missing(nasymp)) nasymp <- object@nasymp
    if (missing(seed)) seed <- object@seed
    if (missing(transform)) transform <- object@transform
    if (missing(transform.data)) transform.data <- object@transform.data
    if (missing(nrbf)) nrbf <- object@nrbf
    if (missing(method)) method <- object@method
    if (missing(lql.frac)) lql.frac <- object@lql.frac
    if (missing(se.par.frac)) se.par.frac <- object@se.par.frac

    nlf(as(object,"pomp"), start=start, est=est, lags=lags, period=period,
      tensor=tensor, nconverge=nconverge, seed=seed, transform=transform,
      transform.data=transform.data, nrbf=nrbf, method=method,
      lql.frac=lql.frac, se.par.frac=se.par.frac, ...)
  }
)

setMethod(
  "logLik",
  signature=signature(object="nlfd_pomp"),
  definition = function(object, ...) {
    object@logql
  }
)

nlf.internal <- function (object, start, est, lags, period, tensor,
  nconverge, nasymp, seed, transform,
  nrbf, method, skip.se, verbose,
  bootstrap, bootsamp, lql.frac, se.par.frac,
  eval.only, transform.data, ...)
{

  ep <- paste0("in ",sQuote("nlf"),": ")
  pompLoad(object,verbose=verbose)

  if (eval.only) est <- character(0)
  if (transform)
    params <- partrans(object,start,dir="toEstimationScale")
  else
    params <- start

  par.index <- which(names(params)%in%est)
  if (length(est)==0) par.index <- integer(0)
  guess <- params[par.index]

  if ((lql.frac<=0)||(lql.frac>=1))
    stop(ep,sQuote("lql.frac")," must be in (0,1)",call.=FALSE)

  if ((se.par.frac<=0)||(se.par.frac>=1))
    stop(ep,sQuote("se.par.frac")," must be in (0,1)",call.=FALSE)

  dt.tol <- 1e-3
  times <- time(object,t0=FALSE)
  t0 <- timezero(object)
  dt <- diff(times)
  if (diff(range(dt))>dt.tol*mean(dt))
    stop(ep,sQuote("nlf")," requires evenly spaced sampling times",call.=FALSE)
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
      transform=transform,
      times=times,
      t0=t0,
      lags=lags,
      period=period,
      tensor=tensor,
      seed=seed,
      transform.data=transform.data,
      nrbf=nrbf,
      verbose=verbose,
      bootstrap=bootstrap,
      bootsamp=bootsamp
    )
    opt <- list(params=params,value=val)
  } else {
    if (method == 'subplex') {
      opt <- subplex::subplex(
        par=guess,
        fn=nlf.objfun,
        object=object,
        params=params,
        par.index=par.index,
        transform=transform,
        times=times,
        t0=t0,
        lags=lags,
        period=period,
        tensor=tensor,
        seed=seed,
        transform.data=transform.data,
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
        gr=NULL,
        method=method,
        object=object,
        params=params,
        par.index=par.index,
        transform=transform,
        times=times,
        t0=t0,
        lags=lags,
        period=period,
        tensor=tensor,
        seed=seed,
        transform.data=transform.data,
        nrbf=nrbf,
        verbose=verbose,
        bootstrap=bootstrap,
        bootsamp=bootsamp,
        control=list(...)
      )
    }

    params[par.index] <- opt$par
    opt$params <- if (transform) partrans(object,params,dir="fromEstimationScale") else params

  }

  opt$Qhat <- matrix(NA,0,0)
  opt$se <- numeric(0)

  ## compute estimated Variance-Covariance matrix of fitted parameters
  fitted <- params[par.index]
  nfitted <- length(fitted)

  if (!skip.se && nfitted>0) {
    Jhat <- matrix(0,nfitted,nfitted)
    Ihat <- Jhat
    f0 <- nlf.lql(
      fitted,
      object=object,
      params=params,
      par.index=par.index,
      transform=transform,
      times=times,
      t0=t0,
      lags=lags,
      period=period,
      tensor=tensor,
      seed=seed,
      transform.data=transform.data,
      nrbf=4,
      verbose=verbose
    )
    F0 <- mean(f0,na.rm=TRUE)

    npts <- length(f0)
    ## Number of lags to use in Newey-West covariance estimator:
    nlags <- round(5*npts^0.2)

    ## find a good epsilon
    h <- se.par.frac
    if (verbose) cat("h in",sQuote("nlf"),"=",h,"\n")
    eps <- rep(h,nfitted)

    for (i in seq_len(nfitted)) {
      Fvals <- rep(0,5)
      Fvals[3] <- F0
      guess <- fitted
      guess[i] <- fitted[i]-sqrt(2)*h*abs(fitted[i])
      Fvals[1] <- mean(
        nlf.lql(
          guess,object=object, params=params, par.index=par.index,
          transform=transform,
          times=times, t0=t0, lags=lags, period=period, tensor=tensor,
          seed=seed, transform.data=transform.data,nrbf=4,
          verbose=verbose
        ),
        na.rm=TRUE
      )
      guess <- fitted
      guess[i] <- fitted[i]-h*abs(fitted[i])
      Fvals[2] <- mean(
        nlf.lql(
          guess,object=object, params=params, par.index=par.index,
          transform=transform,
          times=times, t0=t0, lags=lags, period=period, tensor=tensor,
          seed=seed, transform.data=transform.data, nrbf=4,
          verbose=verbose
        ),
        na.rm=TRUE
      )
      guess <- fitted
      guess[i] <- fitted[i]+h*abs(fitted[i])
      Fvals[4] <- mean(
        nlf.lql(
          guess,object=object, params=params, par.index=par.index,
          transform=transform,
          times=times, t0=t0, lags=lags, period=period, tensor=tensor,
          seed=seed, transform.data=transform.data, nrbf=4,
          verbose=verbose
        ),
        na.rm=TRUE
      )
      guess <- fitted
      guess[i] <- fitted[i]+sqrt(2)*h*abs(fitted[i])
      Fvals[5] <- mean(
        nlf.lql(
          guess,object=object, params=params, par.index=par.index,
          transform=transform,
          times=times, t0=t0, lags=lags, period=period, tensor=tensor,
          seed=seed, transform.data=transform.data, nrbf=4,
          verbose=verbose
        ),
        na.rm=TRUE
      )
      FAILED <- -999999
      Fvals[Fvals < FAILED+10] <- NA
      xvals <- cbind(1,(c(sqrt(2),1,0,1,sqrt(2))*h*fitted[i])^2)
      c2 <- .lm.fit(xvals,Fvals)$coefficients[2]
      eps[i] <- sqrt(abs(lql.frac/c2))
    }

    if (verbose) cat("epsilon in",sQuote("nlf"),"=",t(eps),"\n")

    Imat <- matrix(0,npts,nfitted)
    for (i in seq_len(nfitted)) {
      guess.up <- fitted
      guess.up[i] <- guess.up[i]+eps[i]
      f.up <- nlf.lql(
        guess.up,object=object, params=params, par.index=par.index,
        transform=transform,
        times=times, t0=t0, lags=lags, period=period, tensor=tensor,
        seed=seed, transform.data=transform.data, nrbf=4,
        verbose=verbose
      )
      F.up <- mean(f.up,na.rm=TRUE)

      if (verbose)
        cat("fitted param",i,F.up,"up in",sQuote("nlf"),"\n")

      guess.down <- fitted
      guess.down[i] <- guess.down[i]-eps[i]
      f.down <- nlf.lql(
        guess.down,object=object, params=params, par.index=par.index,
        transform=transform,
        times=times, t0=t0, lags=lags, period=period, tensor=tensor,
        seed=seed, transform.data=transform.data, nrbf=4,
        verbose=verbose
      )
      F.down <- mean(f.down,na.rm=TRUE)

      if (verbose)
        cat("fitted param",i,F.down,"down in",sQuote("nlf"),"\n")

      Jhat[i,i] <- (F.up + F.down-2*F0)/(eps[i]*eps[i])
      Imat[,i] <- (f.up-f.down)/(2*eps[i])
      Ihat[i,i] <- Newey.West(Imat[,i],Imat[,i],nlags)
    }

    for (i in seq_len(nfitted-1)) {
      for (j in seq(from=i+1,to=nfitted,by=1)) {
        guess.uu <- fitted
        guess.uu[i] <- guess.uu[i]+eps[i]
        guess.uu[j] <- guess.uu[j]+eps[j]
        F.uu <- mean(
          nlf.lql(
            guess.uu,object=object, params=params, par.index=par.index,
            transform=transform,
            times=times, t0=t0, lags=lags, period=period, tensor=tensor,
            seed=seed, transform.data=transform.data, nrbf=4,
            verbose=verbose
          ),
          na.rm=TRUE
        )

        guess.ud <- fitted
        guess.ud[i] <- guess.ud[i]+eps[i]
        guess.ud[j] <- guess.ud[j]-eps[j]
        F.ud <- mean(
          nlf.lql(
            guess.ud,object=object, params=params, par.index=par.index,
            transform=transform,
            times=times, t0=t0, lags=lags, period=period, tensor=tensor,
            seed=seed, transform.data=transform.data, nrbf=4,
            verbose=verbose
          ),
          na.rm=TRUE
        )

        guess.du <- fitted
        guess.du[i] <- guess.du[i]-eps[i]
        guess.du[j] <- guess.du[j]+eps[j]
        F.du <- mean(
          nlf.lql(
            guess.du,object=object, params=params, par.index=par.index,
            transform=transform,
            times=times, t0=t0, lags=lags, period=period, tensor=tensor,
            seed=seed, transform.data=transform.data, nrbf=4,
            verbose=verbose
          ),
          na.rm=TRUE
        )

        guess.dd <- fitted
        guess.dd[i] <- guess.dd[i]-eps[i]
        guess.dd[j] <- guess.dd[j]-eps[j]
        F.dd <- mean(
          nlf.lql(
            guess.dd,object=object, params=params, par.index=par.index,
            transform=transform,
            times=times, t0=t0, lags=lags, period=period, tensor=tensor,
            seed=seed, transform.data=transform.data, nrbf=4,
            verbose=verbose
          ),
          na.rm=TRUE
        )

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
    opt$se <- setNames(sqrt(diag(Qhat))/sqrt(npts),names(params)[par.index])
    opt$npts <- npts
  }

  pompUnload(object,verbose=verbose)

  new(
    "nlfd_pomp",
    object,
    params=opt$params,
    transform=transform,
    transform.data=transform.data,
    est=est,
    lags=lags,
    nconverge=nconverge,
    nasymp=nasymp,
    seed=seed,
    period=period,
    tensor=tensor,
    nrbf=nrbf,
    method=method,
    lql.frac=lql.frac,
    se.par.frac=se.par.frac,
    Qhat=opt$Qhat,
    se=opt$se,
    logql=-opt$value
  )
}
