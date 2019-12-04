##' Nonlinear forecasting
##'
##' Parameter estimation by maximum simulated quasi-likelihood.
##'
##' Nonlinear forecasting (NLF) is an \sQuote{indirect inference} method.
##' The NLF approximation to the log likelihood of the data series is computed by simulating data from a model, fitting a nonlinear autoregressive model to the simulated time series, and quantifying the ability of the resulting fitted model to predict the data time series.
##' The nonlinear autoregressive model is implemented as a generalized additive model (GAM), conditional on lagged values, for each observation variable.
##' The errors are assumed multivariate normal.
##'
##' The NLF objective function constructed by \code{nlf_objfun} simulates long time series (\code{nasymp} is the number of observations in the simulated times series), perhaps after allowing for a transient period (\code{ntransient} steps).
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
##' @aliases nlf nlf_objfun nlf_objfun,missing-method nlf_objfun,ANY-method
##'
##' @author Stephen P. Ellner, Bruce E. Kendall, Aaron A. King
##'
##' @family pomp parameter estimation methods
##'
##' @return
##' \code{nlf_objfun} constructs a stateful objective function for NLF estimation.
##' Specfically, \code{nlf_objfun} returns an object of class \sQuote{nlf_objfun}, which is a function suitable for use in an \code{\link{optim}}-like optimizer.
##' In particular, this function takes a single numeric-vector argument that is assumed to contain the parameters named in \code{est}, in that order.
##' When called, it will return the negative log quasilikelihood.
##' It is a stateful function:
##' Each time it is called, it will remember the values of the parameters and its estimate of the log quasilikelihood.
##'
##' @references
##'
##' \Ellner1998
##'
##' \Kendall1999
##'
##' \Kendall2005
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
##' @param nrbf integer scalar;
##' the number of radial basis functions to be used at each lag.
##'
##' @param ti,tf required numeric values.
##' NLF works by generating simulating long time series from the model.
##' The simulated time series will be from \code{ti} to \code{tf}, with the same sampling frequency as the data.
##' \code{ti} should be chosen large enough so that transient dynamics have died away.
##' \code{tf} should be chosen large enough so that sufficiently many data points are available to estimate the nonlinear forecasting model well.
##' An error will be generated unless the data-to-parameter ratio exceeds 10 and
##' a warning will be given if the ratio is smaller than 30.
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
##' @param transform.data optional function.
##' If specified, forecasting is performed using data and model simulations transformed by this function.
##' By default, \code{transform.data} is the identity function,
##' i.e., no transformation is performed.
##' The main purpose of \code{transform.data} is to achieve approximately multivariate normal forecasting errors.
##' If the data are univariate, \code{transform.data} should take a scalar and return a
##' scalar.
##' If the data are multivariate, \code{transform.data} should assume a vector input and return a vector of the same length.
##'
##' @example examples/nlf.R
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
    ti="numeric",
    tf="numeric",
    est="character"
  )
)

setGeneric(
  "nlf_objfun",
  function (data, ...)
    standardGeneric("nlf_objfun")
)

setMethod(
  "nlf_objfun",
  signature=signature(data="missing"),
  definition=function (...) {
    reqd_arg("nlf_objfun","data")
  }
)

setMethod(
  "nlf_objfun",
  signature=signature(data="ANY"),
  definition=function (data, ...) {
    undef_method("nlf_objfun",data)
  }
)

##' @name nlf_objfun-data.frame
##' @aliases nlf_objfun,data.frame-method
##' @rdname nlf
setMethod(
  "nlf_objfun",
  signature=signature(data="data.frame"),
  definition=function (data,
    est = character(0), lags, nrbf = 4, ti, tf,
    seed = NULL, transform.data = identity,
    period = NA, tensor = TRUE, fail.value = NA_real_,
    params, rinit, rprocess, rmeasure,
    ..., verbose = getOption("verbose")) {

    tryCatch(
      nlfof.internal(
        object=data,
        est=est,
        lags=lags,
        nrbf=nrbf,
        ti=ti,
        tf=tf,
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
      error = function (e) pStop("nlf_objfun",conditionMessage(e))
    )

  }
)

##' @name nlf_objfun-pomp
##' @aliases nlf_objfun,pomp-method
##' @rdname nlf
setMethod(
  "nlf_objfun",
  signature=signature(data="pomp"),
  definition=function (data,
    est = character(0), lags, nrbf = 4, ti, tf,
    seed = NULL, transform.data = identity,
    period = NA, tensor = TRUE, fail.value = NA,
    ..., verbose = getOption("verbose")) {

    tryCatch(
      nlfof.internal(
        object=data,
        est=est,
        lags=lags,
        nrbf=nrbf,
        ti=ti,
        tf=tf,
        seed=seed,
        period=period,
        tensor=tensor,
        transform.data=transform.data,
        fail.value=fail.value,
        ...,
        verbose=verbose
      ),
      error = function (e) pStop("nlf_objfun",conditionMessage(e))
    )

  }
)

##' @name nlf_objfun-nlf_objfun
##' @aliases nlf_objfun,nlf_objfun-method
##' @rdname nlf
##' @export
setMethod(
  "nlf_objfun",
  signature=signature(data="nlf_objfun"),
  definition=function (data,
    est, lags, nrbf, ti, tf, seed = NULL,
    period, tensor, transform.data, fail.value,
    ..., verbose = getOption("verbose", FALSE)) {

    if (missing(est)) est <- data@est
    if (missing(lags)) lags <- data@env$lags
    if (missing(nrbf)) nrbf <- data@env$nrbf
    if (missing(ti)) ti <- data@ti
    if (missing(tf)) tf <- data@tf
    if (missing(period)) period <- data@env$period
    if (missing(tensor)) tensor <- data@env$tensor
    if (missing(transform.data)) transform.data <- data@env$transform.data
    if (missing(fail.value)) fail.value <- data@env$fail.value

    nlf_objfun(
      data@env$object,
      est=est,
      lags=lags,
      nrbf=nrbf,
      ti=ti,
      tf=tf,
      seed=seed,
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
  est, lags, nrbf, ti, tf, seed, period, tensor,
  transform.data, fail.value,
  ..., verbose)
{

  verbose <- as.logical(verbose)

  object <- pomp(object, ..., verbose=verbose)

  if (!undefined(object@covar))
    pStop_("NLF is incompatible with time-varying covariates.")

  est <- as.character(est)
  est <- est[nzchar(est)]

  if (!is.numeric(lags) || !all(is.finite(lags)) || any(lags < 1))
    pStop_(sQuote("lags")," must be positive integers.")
  lags <- as.integer(lags)

  if (!is.numeric(nrbf) || length(nrbf) != 1 || !is.finite(nrbf) || nrbf < 4)
    pStop_(sQuote("nrbf")," must be at least 4.")
  nrbf <- as.integer(nrbf)

  if (!(is.logical(period) || is.numeric(period)) ||
      length(period) != 1 || is.infinite(period))
    pStop_(sQuote("period")," must be single finite number or NA.")
  period <- as.numeric(period)
  if (!is.na(period) && period <= 0) period <- NA_real_

  tensor <- as.logical(tensor)

  ## check that times are equally spaced
  dt.tol <- 1e-3
  dts <- diff(time(object,t0=FALSE))
  dt <- mean(dts)
  if (diff(range(dts)) > dt.tol*dt)
    pStop_("NLF requires uniform sampling frequency.")

  t0 <- timezero(object)

  if (!is.numeric(ti) || length(ti) != 1 || !is.finite(ti) || ti < t0)
    pStop_(sQuote("ti")," must be a single numeric value larger than ",sQuote("t0"),".")

  if (!is.numeric(tf) || length(tf) != 1 || !is.finite(tf) || tf <= ti)
    pStop_(sQuote("tf")," must be a single numeric value larger than ",sQuote("ti"),".")

  times <- seq(from=ti,to=tf,by=dt)

  dof <- length(lags)*nrbf
  if (isTRUE(period > 0)) {
    if (tensor) {
      dof <- dof*nrbf
    } else {
      dof <- dof+nrbf
    }
  }

  if (length(times) < 10*dof)
    pStop_("insufficiently long simulated time series: increase ",sQuote("tf"),
      " to at least ",ti+10*dof*dt,".")
  if (length(times) < 30*dof)
    pWarn("nlf_objfun","insufficiently long simulated time series: ",
      "consider increasing ",sQuote("tf")," to ",ti+30*dof*dt," or larger.")

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

  logql <- nlf.lql(object,times=times,lags=lags,nrbf=nrbf,fail.value=fail.value,
    period=period,tensor=tensor,seed=seed,transform.data=transform.data,
    verbose=verbose)

  ofun <- function (par) {
    params[idx] <- par
    coef(object,transform=TRUE) <<- params
    logql <<- nlf.lql(object,times=times,lags=lags,nrbf=nrbf,fail.value=fail.value,
      period=period,tensor=tensor,seed=seed,transform.data=transform.data,
      verbose=verbose)
    if (is.finite(logql) || is.na(fail.value)) -logql else fail.value
  }

  environment(ofun) <- list2env(
    list(object=object,times=times,lags=lags,nrbf=nrbf,fail.value=fail.value,
      period=period,tensor=tensor,seed=seed,transform.data=transform.data,
      params=params,idx=idx,logql=logql,verbose=verbose),
    parent=parent.frame(2)
  )

  new("nlf_objfun",ofun,env=environment(ofun),ti=ti,tf=tf,est=est)

}

## Evaluates the NLF objective function given a POMP object.
##
## Computes the vector of log quasi-likelihood values at the observations
## Note that the log QL itself is returned, not the negative log QL,
## so a large NEGATIVE value is used to flag bad parameters

nlf.lql <- function (object, times, lags, nrbf, period, tensor,
  seed, transform.data, fail.value, verbose)
{

  y <- simulate(object,times=times,seed=seed,format="arrays",verbose=verbose)$obs

  ## Test whether the model time series is valid
  if (!all(is.finite(y))) return(-fail.value)

  dat.mat <- obs(object)
  mod.mat <- array(dim=dim(y)[c(1L,3L)],dimnames=list(rownames(dat.mat),NULL))
  tryCatch(
    {
      dat.mat[,] <- apply(dat.mat,2L,transform.data)
      mod.mat[,] <- apply(y[,1,,drop=FALSE],c(2,3),transform.data)
    },
    error = function (e) pStop("transform.data",conditionMessage(e))
  )

  tryCatch(
    nlf.internal(
      dat.mat=dat.mat,
      dat.times=time(object),
      mod.mat=mod.mat,
      mod.times=times,
      lags=lags,
      nrbf=nrbf,
      period=period,
      tensor=tensor,
      fail.value=fail.value,
      verbose=verbose
    ),
    error = function (e) pStop("nlf.lql",conditionMessage(e))
  )

}

## CENTRAL NLF FUNCTION
##
## Version 1.0, 4 December 2007, S.P. Ellner and Bruce E. Kendall
## Version 1.1, 19 June 2008, S.P. Ellner and A.A. King.
## Version 1.2, 3 September 2018, A.A. King

## Peculiarity of the code: No basis functions involving cross terms of observables.

#####################################################################################
## ARGUMENTS:
## dat.mat = matrix of data time series, nobs x ntimes.data
## mod.mat = matrix of model time series, nobs x ntimes.sim
## lags = vector of lag times for forecasting y(t) = f(y(t-lag1),y(t-lag2),....)+error
## nrbf = number of radial basis functions
## period: period=NA means the model is nonseasonal. period=integer>0 is the period of
##         seasonal forcing in 'real time' (the units of mod.times).
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
##        under the forecasting model based on mod.mat. The approximation used
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

nlf.internal <- function (dat.mat, dat.times, mod.mat, mod.times,
  lags, nrbf, period, tensor, fail.value, verbose) {

  seas <- !is.na(period)
  nvar <- nrow(dat.mat)

  knots <- apply(mod.mat,1L,rbf.knots,n=nrbf)

  dat.lags <- make.lags.nlf(x=dat.mat,times=dat.times,lags=lags)
  sim.lags <- make.lags.nlf(x=mod.mat,times=mod.times,lags=lags)

  dat.pred <- array(dim=c(dim(dat.lags$x),nrbf))
  sim.pred <- array(dim=c(dim(sim.lags$x),nrbf))

  for (i in seq_len(nvar)) {
    dat.pred[,i,,] <- rbf.basis(dat.lags$x[,i,],knots=knots[,i])
    sim.pred[,i,,] <- rbf.basis(sim.lags$x[,i,],knots=knots[,i])
  }

  dd <- dim(dat.pred)
  dim(dat.pred) <- c(dd[1],prod(dd[-1]))
  dd <- dim(sim.pred)
  dim(sim.pred) <- c(dd[1],prod(dd[-1]))

  dat.resp <- dat.lags$y
  sim.resp <- sim.lags$y

  if (seas) {

    knots <- seq(from=-0.1,to=1.1*period,length.out=nrbf)
    dat.seas <- rbf.basis(as.matrix(dat.lags$t %% period),knots=knots)
    sim.seas <- rbf.basis(as.matrix(sim.lags$t %% period),knots=knots)
    dd <- dim(dat.seas)
    dim(dat.seas) <- c(dd[1],prod(dd[-1]))
    dd <- dim(sim.seas)
    dim(sim.seas) <- c(dd[1],prod(dd[-1]))

    if (tensor) {

      ## make gam coefficients time-dependent
      dat.pred <- make.tensorbasis(dat.pred,dat.seas)
      sim.pred <- make.tensorbasis(sim.pred,sim.seas)

    } else {

      ## add time-varying intercept
      dat.pred <- cbind(dat.pred,dat.seas)
      sim.pred <- cbind(sim.pred,sim.seas)


    }
  }

  model.lm <- .lm.fit(sim.pred,sim.resp)

  prediction.err <- dat.resp - dat.pred %*% model.lm$coefficients

  LQL <- dmvnorm(prediction.err,sigma=cov(model.lm$residuals),log=TRUE)
  ## NOTE: This could be improved using GLS.

  LQL <- sum(LQL)

  if (verbose) cat(sprintf("logql = %g\n",LQL))

  if (is.finite(LQL)) LQL else -fail.value

}

## Create design array and response matrix out of lagged variables.
## Data series are in the rows of 'x'.
make.lags.nlf <- function (x, times, lags) {

  nvar <- nrow(x)
  m <- length(lags)     ## number of lag variables per observable
  maxlag <- max(lags)
  n <- ncol(x)-maxlag   ## length of lagged time series
  start <- maxlag+1     ## first predicted observation

  X <- array(dim=c(n,nvar,m))
  for (k in seq_len(length(lags)))
    X[,,k] <- t(x[,seq.int(from=start-lags[k],by=1,length.out=n)])

  y <- t(x[,seq.int(from=start,by=1,length.out=n),drop=FALSE])
  t <- times[seq.int(from=start,by=1,length.out=n)]

  list(x=X,y=y,t=t)

}

rbf.knots <- function (X, n, margin = 0.1) {
  r <- range(X)
  seq(
    from=(1+margin)*r[1]-margin*r[2],
    to=(1+margin)*r[2]-margin*r[1],
    length.out=n
  )
}

rbf.basis <- function (X, knots) {
  X1 <- array(dim=c(dim(X),length(knots)))
  X1[] <- as.numeric(X)-rep(knots,each=length(X))
  s <- diff(range(knots))/length(knots)
  fac <- -1/2/s/s
  exp(fac*X1*X1)
}

## GAUSS trimr function:
## trims n1 rows from the start,
## n2 rows from the end of a matrix or vector
# trimr <- function (a, n1, n2) {
#   a[seq.int(from=n1+1,to=NROW(a)-n2,by=1)]
# }

# Newey.West <- function(x, y, maxlag) {
#   w <- 1-seq_len(maxlag)/(maxlag+1)
#   out <- mean(x*y,na.rm=TRUE)
#   for (i in seq_len(maxlag)) {
#     out <- out+
#       w[i]*mean(trimr(x,i,0)*trimr(y,0,i),na.rm=TRUE)+
#       w[i]*mean(trimr(y,i,0)*trimr(x,0,i),na.rm=TRUE)
#   }
#   out
# }

make.tensorbasis <- function(A,B) {
  if (nrow(A) != nrow(B))
    pStop("make.tensorbasis","incompatible matrices.") # nocov
  nA <- ncol(A)
  nB <- ncol(B)
  Tmat <- matrix(0,nrow(A),nA*nB)
  for (j in seq_len(nB)) {
    Tmat[,seq.int(from=(j-1)*nA+1,to=j*nA,by=1)] <- A*B[,j]
  }
  Tmat
}
