##' Ensemble Kalman filters
##'
##' The ensemble Kalman filter and ensemble adjustment Kalman filter.
##'
##' @name kalman
##' @rdname kalman
##' @include pomp_class.R pomp.R workhorses.R
##' @importFrom stats dnorm rnorm
##' @aliases enkf eakf enkf,ANY-method enkf,missing-method
##' eakf,ANY-method eakf,missing-method
##' @author Aaron A. King
##' @family particle_filter_methods
##' @family estimation_methods
##'
##' @inheritSection pomp Note for Windows users
##' 
##' @inheritParams pomp
##' @param Np the number of particles to use.
##' @param h function returning the expected value of the observation given the
##' state.
##' @param C matrix converting state vector into expected value of the
##' observation.
##' @param R matrix; variance of the measurement noise.
##'
##' @return
##' An object of class \sQuote{kalmand_pomp}.
##'
##' @references
##'
##' \Evensen1994
##'
##' \Anderson2001
##'
##' \Evensen2009
NULL

setClass(
  "kalmand_pomp",
  contains="pomp",
  slots=c(
    Np="integer",
    pred.mean="array",
    filter.mean="array",
    forecast="array",
    cond.logLik="numeric",
    loglik="numeric"
  ),
  prototype=prototype(
    Np=0L,
    pred.mean=array(data=numeric(0),dim=c(0,0)),
    filter.mean=array(data=numeric(0),dim=c(0,0)),
    forecast=array(data=numeric(0),dim=c(0,0)),
    cond.logLik=numeric(0),
    loglik=as.double(NA)
  )
)

setGeneric(
  "enkf",
  function (data, ...)
    standardGeneric("enkf")
)

setGeneric(
  "eakf",
  function (data, ...)
    standardGeneric("eakf")
)

setMethod(
  "eakf",
  signature=signature(data="missing"),
  definition=function (...) {
    reqd_arg("eakf","data")
  }
)

setMethod(
  "eakf",
  signature=signature(data="ANY"),
  definition=function (data, ...) {
    undef_method("eakf",data)
  }
)

setMethod(
  "enkf",
  signature=signature(data="missing"),
  definition=function (...) {
    reqd_arg("enkf","data")
  }
)

setMethod(
  "enkf",
  signature=signature(data="ANY"),
  definition=function (data, ...) {
    undef_method("enkf",data)
  }
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

##' @name enkf-data.frame
##' @aliases enkf,data.frame-method
##' @rdname kalman
##' @export
setMethod(
  "enkf",
  signature=signature(data="data.frame"),
  function (data,
    Np, h, R,
    params, rinit, rprocess,
    ..., verbose = getOption("verbose", FALSE)) {

    tryCatch(
      enkf.internal(
        data,
        Np=Np,
        h=h,
        R=R,
        params=params,
        rinit=rinit,
        rprocess=rprocess,
        ...,
        verbose=verbose
      ),
      error = function (e) pStop("enkf",conditionMessage(e))
    )

  }
)

##' @name enkf-pomp
##' @aliases enkf,pomp-method
##' @rdname kalman
##' @export
setMethod(
  "enkf",
  signature=signature(data="pomp"),
  function (data,
    Np, h, R,
    ..., verbose = getOption("verbose", FALSE)) {

    tryCatch(
      enkf.internal(
        data,
        Np=Np,
        h=h,
        R=R,
        ...,
        verbose=verbose
      ),
      error = function (e) pStop("enkf",conditionMessage(e))
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

##' @name eakf-data.frame
##' @aliases eakf,data.frame-method
##' @rdname kalman
##' @export
setMethod(
  "eakf",
  signature=signature(data="data.frame"),
  function (data,
    Np, C, R,
    params, rinit, rprocess,
    ...,verbose = getOption("verbose", FALSE)) {

    tryCatch(
      eakf.internal(
        data,
        Np=Np,
        C=C,
        R=R,
        params=params,
        rinit=rinit,
        rprocess=rprocess,
        ...,
        verbose=verbose
      ),
      error = function (e) pStop("eakf",conditionMessage(e))
    )

  }
)

##' @name eakf-pomp
##' @aliases eakf,pomp-method
##' @rdname kalman
##' @export
setMethod(
  "eakf",
  signature=signature(data="pomp"),
  function (data,
    Np, C, R,
    ..., verbose = getOption("verbose", FALSE)) {

    tryCatch(
      eakf.internal(
        data,
        Np=Np,
        C=C,
        R=R,
        ...,
        verbose=verbose
      ),
      error = function (e) pStop("eakf",conditionMessage(e))
    )

  }
)

enkf.internal <- function (object,
  Np, h, R,
  ..., verbose) {

  verbose <- as.logical(verbose)

  object <- pomp(object,...,verbose=verbose)

  if (undefined(object@rprocess))
    pStop_(paste(sQuote(c("rprocess")),collapse=", ")," is a needed basic component.")

  Np <- as.integer(Np)
  R <- as.matrix(R)
  params <- coef(object)

  t <- time(object)
  tt <- time(object,t0=TRUE)
  ntimes <- length(t)

  y <- obs(object)
  nobs <- nrow(y)

  pompLoad(object,verbose=verbose)
  on.exit(pompUnload(object,verbose=verbose))

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
      pStop_("degenerate ",sQuote("R"),": ",conditionMessage(e))
    }
  )

  for (k in seq_len(ntimes)) {
    ## advance ensemble according to state process
    X[,] <- rprocess(object,x0=X,t0=tt[k],times=tt[k+1],params=params)

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
    cond.logLik=condlogLik,
    loglik=sum(condlogLik))

}

eakf.internal <- function (object,
  Np, C, R,
  ..., verbose) {

  verbose <- as.logical(verbose)

  object <- pomp(object,...,verbose=verbose)

  if (undefined(object@rprocess))
    pStop_(paste(sQuote(c("rprocess")),collapse=", ")," is a needed basic component.")

  Np <- as.integer(Np)
  R <- as.matrix(R)

  params <- coef(object)

  t <- time(object)
  tt <- time(object,t0=TRUE)

  y <- obs(object)

  pompLoad(object,verbose=verbose)
  on.exit(pompUnload(object,verbose=verbose))

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
      pStop_("degenerate ",sQuote("R"),": ",conditionMessage(e))
    }
  )

  for (k in seq_len(ntimes)) {

    ## advance ensemble according to state process
    X[,] <- rprocess(object,x0=X,t0=tt[k],times=tt[k+1],params=params)

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
    cond.logLik=condlogLik,
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

  list(
    filterMeans=filterMeans,
    predMeans=predMeans,
    forecast=forecast,
    cond.logLik=condlogLik,
    loglik=sum(condlogLik)
  )

}
