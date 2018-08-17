##' Ensemble Kalman filters
##'
##' The ensemble Kalman filter and ensemble adjustment Kalman filter.
##'
##' @name Ensemble Kalman filter methods
##' @rdname kalman
##' @include pomp_class.R
##' @aliases enkf eakf enkf,ANY-method enkf,missing-method
##' eakf,ANY-method eakf,missing-method
##'
##' @param object An object of class \sQuote{pomp} or extending class
##' \sQuote{pomp}.
##' @param params optional named numeric vector containing the parameters at
##' which the filtering should be performed.  By default, \code{params =
##' coef(object)}.
##' @param Np the number of particles to use.
##' @param verbose logical; if \code{TRUE}, progress information is reported.
##' @param h function returning the expected value of the observation given the
##' state.
##' @param C matrix converting state vector into expected value of the
##' observation.
##' @param R matrix; variance of the measurement noise.
##' @param \dots additional arguments, which can be used to modify components of the model.
##'
##' @return An object of class \sQuote{kalmand_pomp}.
##'
##' @author Aaron A. King
##'
##' @seealso \code{\link{pfilter}}, and the tutorials on
##' the \href{https://kingaa.github.io/pomp}{package website}.
##'
##' @references
##' Evensen, G. (1994) Sequential data assimilation with a
##' nonlinear quasi-geostrophic model using Monte Carlo methods to forecast
##' error statistics Journal of Geophysical Research: Oceans 99:10143--10162
##'
##' Evensen, G. (2009) Data assimilation: the ensemble Kalman filter
##' Springer-Verlag.
##'
##' Anderson, J. L. (2001) An Ensemble Adjustment Kalman Filter for Data
##' Assimilation Monthly Weather Review 129:2884--2903
NULL

setClass(
  "kalmand_pomp",
  contains="pomp",
  slots=c(
    Np="integer",
    pred.mean="array",
    filter.mean="array",
    forecast="array",
    cond.loglik="numeric",
    loglik="numeric"
  ),
  prototype=prototype(
    Np=0L,
    pred.mean=array(data=numeric(0),dim=c(0,0)),
    filter.mean=array(data=numeric(0),dim=c(0,0)),
    forecast=array(data=numeric(0),dim=c(0,0)),
    cond.loglik=numeric(0),
    loglik=as.double(NA)
  )
)

setGeneric(
  "enkf",
  function (object, ...)
    standardGeneric("enkf")
)

setGeneric(
  "eakf",
  function (object, ...)
    standardGeneric("eakf")
)

## ENSEMBLE KALMAN FILTER (ENKF)

## Ensemble: $X_t\in \mathbb{R}^{m\times q}$
## Prediction mean: $M_t=\langle X \rangle$
## Prediction variance: $V_t=\langle\langle X \rangle\rangle$
## Forecast: $Y_t=h(X_t)$
## Forecast mean: $N_t=\langle Y \rangle$.
## Forecast variance: $S_t=\langle\langle Y \rangle\rangle$
## State/forecast covariance: $W_t=\langle\langle X,Y\rangle\rangle$
## Kalman gain: $K_t = W_t\,S_t^{-1}$
## New observation: $y_t\in \mathbb{R}^{n\times 1}$
## Updated ensemble: $X^u_{t}=X_t + K_t\,(O_t - Y_t)$
## Filter mean: $m_t=\langle X^u_t \rangle = \frac{1}{q} \sum\limits_{i=1}^q x^{u_i}_t$

##' @name enkf-pomp
##' @aliases enkf,pomp-method
##' @rdname kalman
setMethod(
  "enkf",
  signature=signature(object="pomp"),
  function (object, params, Np, h, R, ...,
    verbose = getOption("verbose", FALSE)) {

    enkf.internal(
      object=object,
      params=params,
      h=h,
      R=R,
      Np=Np,
      ...,
      verbose=verbose
    )

  }
)

## ENSEMBLE ADJUSTMENT KALMAN FILTER (EAKF)

## Ensemble: $X_t\in \mathbb{R}^{m\times q}$
## Prediction mean: $M_t=\langle X \rangle$ (ensemble average).
## Prediction variance: $V_t=\langle\langle X \rangle\rangle$ (ensemble variance).
## SVD of prediction variance: $V = Q_{V}\,D_{V}\,Q_{V}^T$
## Another SVD: $U=D_V^{1/2}\,Q_V^T\,C^T\,R^{-1}\,C\,Q_V\,D_V^{1/2}=Q_U\,D_U\,Q_U^T$
## Adjustment: $B=Q_V\,D_V^{1/2}\,Q_U\,(I+D_U)^{-1/2}\,D_V^{-1/2}\,Q_V^T$
## Kalman gain: $K=B\,V\,B^T\,C^T\,R^{-1}$
## Filter mean: $m_t=M_t+K_t\,(y_t-C\,M_t)$
## Updated ensemble: $x_{t}=B\,(X_t-M_t\,\mathbb{1})+m_t\,\mathbb{1}$

##' @name eakf-pomp
##' @aliases eakf,pomp-method
##' @rdname kalman
setMethod(
  "eakf",
  signature=signature(object="pomp"),
  function (object, params, Np, C, R, ...,
    verbose = getOption("verbose", FALSE)) {

    eakf.internal(
      object=object,
      params=params,
      C=C,
      R=R,
      Np=Np,
      ...,
      verbose=verbose
    )
  }
)

setMethod(
  "eakf",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("eakf"),": ",sQuote("object")," is a required argument",call.=FALSE)
  }
)

setMethod(
  "eakf",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("eakf")," is not defined when ",sQuote("object")," is of class ",sQuote(class(object)),call.=FALSE)
  }
)

setMethod(
  "enkf",
  signature=signature(object="missing"),
  definition=function (...) {
    stop("in ",sQuote("enkf"),": ",sQuote("object")," is a required argument",call.=FALSE)
  }
)

setMethod(
  "enkf",
  signature=signature(object="ANY"),
  definition=function (object, ...) {
    stop(sQuote("enkf")," is not defined when ",sQuote("object")," is of class ",sQuote(class(object)),call.=FALSE)
  }
)

enkf.internal <- function (object, params, h, R, Np, ..., verbose) {

  ep <- paste0("in ",sQuote("enkf"),": ")

  Np <- as.integer(Np)
  R <- as.matrix(R)
  verbose <- as.logical(verbose)

  object <- pomp(object,...)

  if (missing(params)) params <- coef(object)
  if (is.list(params)) params <- unlist(params)

  t <- time(object)
  tt <- time(object,t0=TRUE)
  ntimes <- length(t)

  y <- obs(object)
  nobs <- nrow(y)

  Y <- array(dim=c(nobs,Np),dimnames=list(variable=rownames(y),rep=NULL))
  X <- rinit(object,params=params,nsim=Np)
  nvar <- nrow(X)

  filterMeans <- array(dim=c(nvar,ntimes),dimnames=list(variable=rownames(X),time=t))
  predMeans <- array(dim=c(nvar,ntimes),dimnames=list(variable=rownames(X),time=t))
  forecast <- array(dim=c(nobs,ntimes),dimnames=dimnames(y))
  condlogLik <- numeric(ntimes)

  sqrtR <- tryCatch(
    t(chol(R)),                     # t(sqrtR)%*%sqrtR == R
    error = function (e) {
      stop(ep,"degenerate ",sQuote("R"),": ",conditionMessage(e),call.=FALSE)
    }
  )

  for (k in seq_len(ntimes)) {
    ## advance ensemble according to state process
    X[,] <- rprocess(object,params=params,xstart=X,times=tt[c(k,k+1)],offset=1)

    predMeans[,k] <- pm <- rowMeans(X) # prediction mean
    Y[,] <- apply(X,2,h)               # ensemble of forecasts
    ym <- rowMeans(Y)                  # forecast mean

    X <- X-pm
    Y <- Y-ym

    fv <- tcrossprod(Y)/(Np-1)+R    # forecast variance
    vyx <- tcrossprod(Y,X)/(Np-1)   # forecast/state covariance

    svdS <- svd(fv,nv=0)            # singular value decomposition
    Kt <- svdS$u%*%(crossprod(svdS$u,vyx)/svdS$d) # transpose of Kalman gain
    Ek <- sqrtR%*%matrix(rnorm(n=nobs*Np),nobs,Np) # artificial noise
    resid <- y[,k]-ym

    X <- X+pm+crossprod(Kt,resid-Y+Ek)

    condlogLik[k] <- sum(dnorm(x=crossprod(svdS$u,resid),mean=0,sd=sqrt(svdS$d),log=TRUE))
    filterMeans[,k] <- rowMeans(X)  # filter mean
    forecast[,k] <- ym

  }

  new("kalmand_pomp",object,Np=Np,
    filter.mean=filterMeans,
    pred.mean=predMeans,
    forecast=forecast,
    cond.loglik=condlogLik,
    loglik=sum(condlogLik))
}

eakf.internal <- function (object, params, C, R, Np, ..., verbose) {

  ep <- paste0("in ",sQuote("eakf"),": ")

  Np <- as.integer(Np)
  R <- as.matrix(R)
  verbose <- as.logical(verbose)

  object <- pomp(object,...)

  if (missing(params)) params <- coef(object)
  if (is.list(params)) params <- unlist(params)

  t <- time(object)
  tt <- time(object,t0=TRUE)

  y <- obs(object)

  X <- rinit(object,params=params,nsim=Np)

  nvar <- nrow(X)
  ntimes <- length(t)
  nobs <- nrow(y)

  filterMeans <- array(dim=c(nvar,ntimes),
    dimnames=list(variable=rownames(X),time=t))
  predMeans <- array(dim=c(nvar,ntimes),
    dimnames=list(variable=rownames(X),time=t))
  forecast <- array(dim=c(nobs,ntimes),dimnames=dimnames(y))
  condlogLik <- numeric(ntimes)

  ri <- tryCatch(
    solve(R),
    error = function (e) {
      stop(ep,"degenerate ",sQuote("R"),": ",conditionMessage(e),call.=FALSE)
    }
  )

  for (k in seq_len(ntimes)) {

    ## advance ensemble according to state process
    X[,] <- rprocess(object,xstart=X,params=params,times=tt[c(k,k+1)],offset=1)

    predMeans[,k] <- pm <- rowMeans(X) # prediction mean
    X <- X-pm

    pv <- tcrossprod(X)/(Np-1)      # prediction variance
    svdV <- svd(pv,nv=0)

    forecast[,k] <- C %*% pm        # forecast (observables)
    resid <- y[,k]-forecast[,k]     # forecast error

    ww <- t(C%*%svdV$u)
    w <- crossprod(ww,svdV$d*ww)+R  # forecast variance
    svdW <- svd(w,nv=0)

    condlogLik[k] <- sum(dnorm(x=crossprod(svdW$u,resid),mean=0,
      sd=sqrt(svdW$d),log=TRUE))

    u <- sqrt(svdV$d)*ww
    u <- tcrossprod(u%*%ri,u)
    svdU <- svd(u,nv=0)

    ## adjustment
    b <- svdV$u%*%(sqrt(svdV$d)*svdU$u)%*%
      (1/sqrt(1+svdU$d)/sqrt(svdV$d)*t(svdV$u))

    K <- tcrossprod(b%*%pv,b)%*%crossprod(C,ri)   # Kalman gain

    filterMeans[,k] <- fm <- pm+K%*%resid         # filter mean

    X[,] <- b%*%X+fm[,]
  }

  new("kalmand_pomp",object,Np=Np,
    filter.mean=filterMeans,
    pred.mean=predMeans,
    forecast=forecast,
    cond.loglik=condlogLik,
    loglik=sum(condlogLik))
}

## Basic Kalman filter. Currently for internal consumption only.
## X(t) ~ MVN(A*X(t-1),Q)
## Y(t) ~ MVN(C*X(t),R)
kalmanFilter <- function (t, y, X0, A, Q, C, R) {
  N <- ncol(y)
  nvar <- length(X0)
  nobs <- nrow(y)
  filterMeans <- array(dim=c(nvar,N),
    dimnames=list(variable=names(X0),time=NULL))
  predMeans <- filterMeans
  forecast <- array(dim=c(nobs,N),dimnames=dimnames(y))
  condlogLik <- numeric(N)
  ri <- solve(R)
  cric <- crossprod(C,ri)%*%C
  fm <- X0
  fv <- matrix(0,nvar,nvar)
  for (k in seq_along(t)) {
    predMeans[,k] <- pm <- A%*%fm      # prediction mean
    pv <- A%*%tcrossprod(fv,A)+Q       # prediction variance
    svdV <- svd(pv,nv=0)
    resid <- y[,k]-C%*%pm              # forecast error
    w <- tcrossprod(C%*%pv,C)+R        # forecast variance
    svdW <- svd(w,nv=0)
    condlogLik[k] <- sum(dnorm(x=crossprod(svdW$u,resid),mean=0,
      sd=sqrt(svdW$d),log=TRUE))
    pvi <- svdV$u%*%(t(svdV$u)/svdV$d) # prediction precision
    fvi <- pvi+cric                    # filter precision
    svdv <- svd(fvi,nv=0)
    fv <- svdv$u%*%(t(svdv$u)/svdv$d)  # filter variance
    K <- fv%*%crossprod(C,ri)          # Kalman gain
    filterMeans[,k] <- fm <- pm+K%*%resid  # filter mean
    forecast[,k] <- C %*% pm
  }
  list(filterMeans=filterMeans,
    predMeans=predMeans,
    forecast=forecast,
    cond.loglik=condlogLik,
    loglik=sum(condlogLik))
}
