##' Nonlinear forecasting
##'
##' Parameter estimation my maximum simulated quasi-likelihood.
##'
##' NLF is an \sQuote{indirect inference} method.
##' The NLF approximation to the log likelihood of the data series is computed by simulating data from a model, fitting a nonlinear autoregressive model to the simulated time series, and quantifying the ability of the resulting fitted model to predict the data time series.
##' The nonlinear autoregressive model is implemented as a generalized additive model (GAM), conditional on lagged values, for each observation variable.
##' The errors are assumed multivariate normal.
##'
##' The NLF objective function constructed by \code{nlf.objfun} simulates long time series (\code{nasymp} is the number of observations in the simulated times series), perhaps after allowing for a transient period (\code{ntransient} steps).
##' It then fits the GAM for the chosen lags to the simulated time series.
##' Finally, it computes the quasi-likelihood of the data under the fitted GAM.
##'
##' NLF assumes that the observation frequency (equivalently the time between successive observations) is uniform.
##'
##' @name nlf
##' @rdname nlf
##' @include pomp.R simulate.R probe_match.R
##'
##' @importFrom stats .lm.fit optim setNames dnorm .lm.fit sd cov
##' @importFrom mvtnorm dmvnorm
##'
##' @aliases nlf nlf.objfun nlf.objfun,missing-method nlf.objfun,ANY-method
##'
##' @author Stephen P. Ellner, Bruce E. Kendall, Aaron A. King
##'
##' @family \pkg{pomp} parameter estimation methods
##'
##' @return
##' \code{nlf.objfun} constructs a stateful objective function for NLF estimation.
##' Specfically, \code{nlf.objfun} returns an object of class \sQuote{nlf_objfun}, which is a function suitable for use in an \code{\link{optim}}-like optimizer.
##' In particular, this function takes a single numeric-vector argument that is assumed to contain the parameters named in \code{est}, in that order.
##' When called, it will return the negative log quasilikelihood.
##' It is a stateful function:
##' Each time it is called, it will remember the values of the parameters and its estimate of the log quasilikelihood.
##'
##' @references
##' Ellner, S. P., Bailey, B. A., Bobashev, G. V., Gallant, A. R., Grenfell, B. T. and Nychka D. W. (1998)
##' Noise and nonlinearity in measles epidemics:
##' combining mechanistic and statistical approaches to population modeling.
##' \emph{American Naturalist} \bold{151}, 425--440.
##'
##' Kendall, B. E., Briggs, C. J., Murdoch, W. W., Turchin, P., Ellner, S. P., McCauley, E., Nisbet, R. M. and Wood S. N. (1999)
##' Why do populations cycle?
##' A synthesis of statistical and mechanistic modeling approaches.
##' \emph{Ecology} \bold{80}, 1789--1805.
##'
##' Kendall, B. E., Ellner, S. P., McCauley, E., Wood, S. N., Briggs, C. J., Murdoch, W. W. and Turchin, P. (2005)
##' Population cycles in the pine looper moth (\emph{Bupalus piniarius}):
##' dynamical tests of mechanistic hypotheses.
##' \emph{Ecological Monographs} \bold{75}, 259--276.
##'
##' @section Periodically-forced systems (seasonality):
##' Unlike other \pkg{pomp} estimation methods, NLF cannot accommodate general time-dependence in the model via explicit time-dependence or dependence on time-varying covariates.
##' However, NLF can accommodate periodic forcing.
##' It does this by including forcing phase as a predictor in the nonlinear autoregressive model.
##' To accomplish this, one sets \code{period} to the period of the forcing (a positive numerical value).
##' In this case, if \code{tensor = FALSE}, the effect is to add a periodic intercept in the autoregressive model.
##' If \code{tensor = TRUE}, by contrast, the fitted model includes time-periodic coefficients,
##' constructed using tensor products of basis functions of observables with
##' basis functions of time.
##'
##' @inheritSection objfun Important Note
##'
##' @param lags A vector specifying the lags to use when constructing the nonlinear autoregressive prediction model.
##' The first lag is the prediction interval.
##'
##' @param period numeric; \code{period=NA} means the model is nonseasonal.
##' \code{period > 0} is the period of seasonal forcing.
##' \code{period <= 0} is equivalent to \code{period = NA}.
##'
##' @param tensor logical;
##' if FALSE, the fitted model is a generalized additive model with time mod period as one of the predictors,
##' i.e., a gam with time-varying intercept.
##' If TRUE, the fitted model is a gam with lagged state variables as predictors and time-periodic coefficients, constructed using tensor products of basis functions of state variables with basis
##' functions of time.
##'
##' @param ntransient number of timesteps to be discarded from the model simulation.
##'
##' @param nasymp number of asymptotic timesteps to be recorded from the model simulation.
##'
##' @param nrbf integer scalar;
##' the number of radial basis functions to be used at each lag.
##'
##' @param transform.data optional function.
##' If specified, forecasting is performed using data and model simulations transformed by this function.
##' By default, \code{transform.data} is the identity function,
##' i.e., no transformation is performed.
##' The main purpose of \code{transform.data} is to achieve approximately multivariate normal forecasting errors.
##' If the data are univariate, \code{transform.data} should take a scalar and return a
##' scalar.
##' If the data are multivariate, \code{transform.data} should assume a vector input and return a vector of the same length.
##'
##' @inheritParams probe.match
##' @inheritParams pomp
##'
NULL

## nonlinear forecasting
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
##
## v 2.0, 3 September 2018, by Aaron King
## - removes bootstrapping and standard errors
## - incorporates the new stateful-function approach of pomp version 2

setClass(
  "nlf_objfun",
  contains="function",
  slots=c(
    env="environment",
    ntransient="integer",
    nasymp="integer"
  )
)

setGeneric(
  "nlf.objfun",
  function (data, ...)
    standardGeneric("nlf.objfun")
)

setMethod(
  "nlf.objfun",
  signature=signature(data="missing"),
  definition=function (...) {
    reqd_arg("nlf.objfun","data")
  }
)

setMethod(
  "nlf.objfun",
  signature=signature(data="ANY"),
  definition=function (data, ...) {
    undef_method("nlf.objfun",data)
  }
)

##' @name logLik-nlf_objfun
##' @aliases logLik,nlf_objfun-method
##' @rdname loglik
##'
##' @return
##' When \code{object} is an NLF objective function, i.e., the result of a call to \code{nlf.objfun},
##' \code{logLik} retrieves the \dQuote{quasi log likelihood} (see \code{\link{nlf}}).
##'
##' @export
setMethod(
  "logLik",
  signature=signature(object="nlf_objfun"),
  definition = function(object, ...) {
    object@env$loglik
  }
)

##' @name nlf.objfun-data.frame
##' @aliases nlf.objfun,data.frame-method
##' @rdname nlf
setMethod(
  "nlf.objfun",
  signature=signature(data="data.frame"),
  definition=function (data,
    est = character(0), lags, nrbf = 4, ntransient, nasymp,
    seed = NULL, transform.data = identity,
    period = NA, tensor = FALSE, fail.value = NA_real_,
    params, rinit, rprocess, rmeasure,
    ..., verbose = getOption("verbose")) {

    tryCatch(
      nlfof.internal(
        object=data,
        est=est,
        lags=lags,
        nrbf=nrbf,
        ntransient=ntransient,
        nasymp=nasymp,
        seed=seed,
        period=period,
        tensor=tensor,
        transform.data=transform.data,
        fail.value=fail.value,
        params=params,
        rinit=rinit,
        rprocess=rprocess,
        rmeasure=rmeasure,
        ...,
        verbose=verbose
      ),
      error = function (e) pStop("nlf.objfun",conditionMessage(e))
    )

  }
)

##' @name nlf.objfun-pomp
##' @aliases nlf.objfun,pomp-method
##' @rdname nlf
setMethod(
  "nlf.objfun",
  signature=signature(data="pomp"),
  definition=function (data,
    est = character(0), lags, nrbf = 4, ntransient, nasymp,
    seed = NULL, transform.data = identity,
    period = NA, tensor = FALSE, fail.value = NA,
    ..., verbose = getOption("verbose")) {

    tryCatch(
      nlfof.internal(
        object=data,
        est=est,
        lags=lags,
        nrbf=nrbf,
        ntransient=ntransient,
        nasymp=nasymp,
        seed=seed,
        period=period,
        tensor=tensor,
        transform.data=transform.data,
        fail.value=fail.value,
        ...,
        verbose=verbose
      ),
      error = function (e) pStop("nlf.objfun",conditionMessage(e))
    )

  }
)

##' @name nlf.objfun-nlf_objfun
##' @aliases nlf.objfun,nlf_objfun-method
##' @rdname nlf
##' @export
setMethod(
  "nlf.objfun",
  signature=signature(data="nlf_objfun"),
  definition=function (data,
    est, lags, nrbf, ntransient, nasymp, seed = NULL,
    period, tensor, transform.data, fail.value,
    ..., verbose = getOption("verbose", FALSE)) {

    if (missing(est)) est <- data@env$est
    if (missing(lags)) lags <- data@env$lags
    if (missing(nrbf)) nrbf <- data@env$nrbf
    if (missing(ntransient)) ntransient <- data@ntransient
    if (missing(nasymp)) nasymp <- data@nasymp
    if (missing(period)) period <- data@env$period
    if (missing(tensor)) tensor <- data@env$tensor
    if (missing(transform.data)) transform.data <- data@env$transform.data
    if (missing(fail.value)) fail.value <- data@env$fail.value

    nlf.objfun(
      data@env$object,
      est=est,
      lags=lags,
      nrbf=nrbf,
      ntransient=ntransient,
      nasymp=nasymp,
      period=period,
      tensor=tensor,
      transform.data=transform.data,
      fail.value=fail.value,
      ...,
      verbose=verbose
    )

  }
)

nlfof.internal <- function (object,
  est, lags, nrbf, ntransient, nasymp, seed, period, tensor,
  transform.data, fail.value,
  ..., verbose)
{

  verbose <- as.logical(verbose)

  object <- pomp(object, ..., verbose=verbose)

  est <- as.character(est)
  est <- est[nzchar(est)]

  if (!is.numeric(lags) || !all(is.finite(lags)) || any(lags < 1))
    pStop_(sQuote("lags")," must be positive integers.")
  lags <- as.integer(lags)

  if (!is.numeric(nrbf) || length(nrbf) != 1 || !is.finite(nrbf) || nrbf < 4)
    pStop_(sQuote("nrbf")," must be at least 4.")
  nrbf <- as.integer(nrbf)

  if (!(is.logical(period) || is.numeric(period)) || length(period) != 1 || is.infinite(period))
    pStop_(sQuote("period")," must be single finite number or NA.")
  period <- as.numeric(period)
  if (!is.na(period) && period <= 0) period <- NA_real_

  tensor <- as.logical(tensor)

  if (!is.numeric(nasymp) || length(nasymp) != 1 || !is.finite(nasymp) || nasymp < 0)
    pStop_(sQuote("nasymp")," must be a nonnegative integer (and should be large).")
  nasymp <- as.integer(nasymp)

  if (!is.numeric(ntransient) || length(ntransient) != 1 || !is.finite(ntransient) || ntransient < max(lags)+1)
    pStop_(sQuote("ntransient")," must be larger than the maximum lag ","(and should be much larger).")
  ntransient <- as.integer(ntransient)

  transform.data <- match.fun(transform.data)

  if (!(is.logical(fail.value) || is.numeric(fail.value)) || length(fail.value) != 1)
    pStop_(sQuote("fail.value")," should be a single (large) number or ",sQuote("NA"),".")
  fail.value <- as.numeric(fail.value)
  if (isTRUE(fail.value < 1000))
    pWarn("nlf",sQuote("fail.value")," should be a large number or ",sQuote("NA"),".")

  pompLoad(object,verbose=verbose)

  params <- coef(object,transform=TRUE)

  idx <- match(est,names(params))
  if (any(is.na(idx))) {
    missing <- est[is.na(idx)]
    pStop_("parameter",ngettext(length(missing),"","s")," ",
      paste(sQuote(missing),collapse=",")," not found in ",sQuote("params"),".")
  }

  ## check that times are equally spaced
  dt.tol <- 1e-3
  times <- time(object,t0=FALSE)
  dt <- diff(times)
  if (diff(range(dt))>dt.tol*mean(dt))
    pStop_("NLF requires uniform sampling frequency.")

  dt <- times[2]-times[1]

  t0 <- timezero(object)
  times <- seq(from=t0+ntransient*dt,length.out=nasymp,by=dt)

  loglik <- nlf.lql(object,times=times,lags=lags,nrbf=nrbf,fail.value=fail.value,
    period=period,tensor=tensor,seed=seed,transform.data=transform.data)

  ofun <- function (par) {
    params[idx] <- par
    coef(object,transform=TRUE) <<- params
    loglik <<- nlf.lql(object,times=times,lags=lags,nrbf=nrbf,fail.value=fail.value,
      period=period,tensor=tensor,seed=seed,transform.data=transform.data)
    if (is.finite(loglik) || is.na(fail.value)) -loglik else fail.value
  }

  environment(ofun) <- list2env(
    list(object=object,times=times,lags=lags,nrbf=nrbf,fail.value=fail.value,
      period=period,tensor=tensor,seed=seed,transform.data=transform.data,
      params=params,idx=idx,loglik=loglik),
    parent=parent.frame(2)
  )

  new("nlf_objfun",ofun,env=environment(ofun),ntransient=ntransient,nasymp=nasymp)

}

## Evaluates the NLF objective function given a POMP object.
##
## Computes the vector of log quasi-likelihood values at the observations
## Note that the log QL itself is returned, not the negative log QL,
## so a large NEGATIVE value is used to flag bad parameters

nlf.lql <- function (object,
  times, lags, nrbf, period, tensor, seed, transform.data, fail.value)
{

  y <- simulate(object,times=times,seed=seed,format="arrays")$obs

  ## Test whether the model time series is valid
  if (!all(is.finite(y))) return(-fail.value)

  data.mat <- obs(object)
  data.mat[,] <- apply(data.mat,2L,transform.data)

  model.mat <- array(
    dim=dim(y)[c(1L,3L)],
    dimnames=list(rownames(data.mat),NULL)
  )
  model.mat[,] <- apply(y[,1,,drop=FALSE],c(2,3),transform.data)

  tryCatch(
    nlf.internal(
      data.mat=data.mat,
      data.times=time(object),
      model.mat=model.mat,
      model.times=times,
      lags=lags,
      nrbf=nrbf,
      period=period,
      tensor=tensor,
      fail.value=fail.value
    ),
    error = function (e) pStop("nlf.lql",conditionMessage(e))
  )

}

## Key 'nlf' function
nlf.internal <- function (data.mat, data.times, model.mat, model.times,
  lags, nrbf, period, tensor, fail.value) {

  ## Version 1.0, 4 December 2007, S.P. Ellner and Bruce E. Kendall
  ## Version 1.1, 19 June 2008, S.P. Ellner and A.A. King.
  ## Version 1.2, 3 September 2018, A.A. King

  ## Peculiarity of the code: No basis functions involving cross terms of observables.

  #####################################################################################
  ## ARGUMENTS:
  ## data.mat = matrix of data time series, nobs x ntimes.data
  ## model.mat = matrix of model time series, nobs x ntimes.sim
  ## lags = vector of lag times for forecasting y(t) = f(y(t-lag1),y(t-lag2),....)+error
  ## nrbf = number of radial basis functions
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

  FAILED = -fail.value

  nvar <- nrow(data.mat)
  multivar <- (nvar > 1)
  seas <- (!is.na(period) && period > 0)

  if (seas) {
    seas.sim <- model.times %% period
    seas.data <- data.times %% period
  } else {
    seas.sim <- NULL
    seas.data <- NULL
  }

  ## do a lagged embedding for observation variable 1
  data.ts <- data.mat[1L,]
  model.ts <- model.mat[1L,]

  ## Set up the RBF knots
  xm <- diff(range(model.ts))
  rbf.knots <- min(model.ts)+seq(-0.1,1.1,length=nrbf)*xm
  s <- 0.3*xm
  fac <- -1/(2*s^2)
  if (!is.finite(fac) || fac == 0) return(FAILED)

  ## Lag the data and set up the predicted values & seasonal indices
  Lags.model <- make.lags.nlf(model.ts,lags=lags,cov=seas.sim)
  Lags.data <- make.lags.nlf(data.ts,lags=lags,cov=seas.data)

  data.pred <- matrix(Lags.data$y,ncol=1)
  model.pred <- matrix(Lags.model$y,ncol=1)

  rbfbasis.model <- make.rbfbasis(Lags.model$x,knots=rbf.knots,fac=fac)
  rbfbasis.data <- make.rbfbasis(Lags.data$x,knots=rbf.knots,fac=fac)

  if (seas) {
    ## Set up the RBF knots
    rbf.cov.knots <- seq(-0.1,1.1,length=nrbf)*period
    s <- 0.3*period
    fac.cov <- -1/(2*s^2)
    if (!is.finite(fac.cov) || fac == 0) return(FAILED)

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
      Lags.model <- make.lags.nlf(model.ts,lags=lags,cov=seas.sim)
      Lags.data <- make.lags.nlf(data.ts,lags=lags,cov=seas.data)

      data.pred <- cbind(data.pred,Lags.data$y)
      model.pred <- cbind(model.pred,Lags.model$y)
      rbfbasis.model <- cbind(rbfbasis.model,make.rbfbasis(Lags.model$x,knots=rbf.knots,fac=fac))
      rbfbasis.data <- cbind(rbfbasis.data,make.rbfbasis(Lags.data$x,knots=rbf.knots,fac=fac))

    }
  }

  if (seas) {
    if (tensor) {
      ## make gam coefficients time-dependent
      rbfbasis.model <- make.tensorbasis.nlf(rbfbasis.model,rbfbasis.cov.model)
      rbfbasis.data <- make.tensorbasis.nlf(rbfbasis.data,rbfbasis.cov.data)
    } else {
      ## add time-varying intercept
      rbfbasis.model <- cbind(rbfbasis.model,rbfbasis.cov.model)
      rbfbasis.data <- cbind(rbfbasis.data,rbfbasis.cov.data)
    }
  }

  prediction.errors <- matrix(0,nrow(data.pred),nvar)
  model.residuals <- matrix(0,nrow(model.pred),nvar)

  for (jvar in seq_len(nvar)) {
    model.lm <- .lm.fit(rbfbasis.model,model.pred[,jvar])
    model.residuals[,jvar] <- model.lm$residuals
    ck <- model.lm$coefficients
    fitted.data <- rbfbasis.data%*%matrix(ck,ncol=1)
    prediction.errors[,jvar] <- data.pred[,jvar]-fitted.data
  }

  if (multivar) {
    sigma.model <- cov(model.residuals)
    LQL <- dmvnorm(prediction.errors,sigma=sigma.model,log=TRUE)
    ## NOTE: This could be improved using GLS.
  } else {
    sigma.model <- sd(model.residuals[,1])
    LQL <- dnorm(prediction.errors[,1],mean=0,sd=sigma.model,log=TRUE)
  }

  sum(LQL)
}

make.lags.nlf <- function(x, lags, cov = NULL, nobs = 10000) {
  x <- as.matrix(x)
  xd <- ncol(x)
  m <- length(lags)
  N <- min(nobs,nrow(x)-max(lags))
  n <- min(nobs,N)
  if (N > nobs) pWarn("make.lags.nlf","series length truncated to default.")
  start <- max(lags)+1
  temp <- matrix(0,ncol=xd*length(lags),nrow=n)
  for (k in seq_len(length(lags))) {
    a <- start-lags[k]
    b <- a + n - 1
    temp[,(1:xd)+(k-1)*xd] <- x[(a:b),]
  }
  a <- start
  b <- a + n - 1
  if (xd == 1)
    lab <- format(paste0("lag",rep(lags,rep(xd,length(lags)))))
  else
    lab <- format(paste0(rep(seq_len(xd),length(lags)),"lag",rep(lags,rep(xd,length(lags)))))
  dimnames(temp) <- list(NULL,lab)
  skip <- NA
  if (!is.null(cov)) {
    cov <- as.matrix(cov)
    cov <- cov[a:b,,drop=FALSE]
    skip <- seq_len(ncol(cov))+m*xd
  }
  if (xd == 1)
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
    for (j in seq(from=2,to=nknots,by=1)) {
      X1 <- cbind(X1,X-knots[j])
    }
  }
  exp(fac*(X1^2))
}

## GAUSS trimr function:
## trims n1 rows from the start,
## n2 rows from the end of a matrix or vector
trimr <- function (a, n1, n2) {
  a[seq.int(from=n1+1,to=NROW(a)-n2,by=1)]
}

Newey.West <- function(x, y, maxlag) {
  w <- 1-seq_len(maxlag)/(maxlag+1)
  out <- mean(x*y,na.rm=TRUE)
  for (i in seq_len(maxlag)) {
    out <- out+
      w[i]*mean(trimr(x,i,0)*trimr(y,0,i),na.rm=TRUE)+
      w[i]*mean(trimr(y,i,0)*trimr(x,0,i),na.rm=TRUE)
  }
  out
}

make.tensorbasis.nlf <- function(A,B) {
  if(nrow(A)!=nrow(B)) pStop("make.tensorbasis.nlf","incompatible matrices.")
  ncol.A <- ncol(A)
  ncol.B <- ncol(B)
  Tmat <- matrix(0,nrow(A),ncol.A*ncol.B)
  for (i in seq_len(ncol.A)) {
    start <- (i-1)*ncol.B
    for (j in seq_len(ncol.B)) {
      Tmat[,start+j] <- A[,i]*B[,j]
    }
  }
  Tmat
}
