##' Ensemble Kalman filters
##'
##' The ensemble Kalman filter and ensemble adjustment Kalman filter.
##'
##' @name kalman
##' @rdname kalman
##' @include pomp_class.R pomp.R workhorses.R
##' @importFrom stats dnorm rnorm
##' @aliases enkf,ANY-method enkf,missing-method
##' eakf,ANY-method eakf,missing-method
##' @author Aaron A. King
##' @concept Kalman filter
##' @seealso \code{\link{kalmanFilter}}
##' @family particle filter methods
##' @family elementary algorithms
##' @inheritSection pomp Note for Windows users
##' @inheritParams pfilter
##' @inheritParams pomp
##' @param Np integer; the number of particles to use, i.e., the size of the ensemble.
##' 
##' @return
##' An object of class \sQuote{kalmand_pomp}.
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
    R=array(data=numeric(0L),dim=c(0L,0L)),
    pred.mean=array(data=numeric(0L),dim=c(0L,0L)),
    filter.mean=array(data=numeric(0L),dim=c(0L,0L)),
    forecast=array(data=numeric(0L),dim=c(0L,0L)),
    cond.logLik=numeric(0L),
    loglik=as.double(NA)
  )
)

setClass(
  "eakfd_pomp",
  contains="kalmand_pomp",
  slots=c(
    Cmatrix="array"
  ),
  prototype=prototype(
    Cmatrix=array(data=numeric(0L),dim=c(0L,0L))
  )
)

setGeneric(
  "enkf",
  function (data, ...)
    standardGeneric("enkf")
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

##' @rdname kalman
##' @aliases enkf
##' @export
setMethod(
  "enkf",
  signature=signature(data="data.frame"),
  function (data, Np,
    params, rinit, rprocess, emeasure, vmeasure,
    ..., verbose = getOption("verbose", FALSE)) {

    tryCatch(
      enkf_internal(
        data,
        Np=Np,
        params=params,
        rinit=rinit,
        rprocess=rprocess,
        emeasure=emeasure,
        vmeasure=vmeasure,
        ...,
        verbose=verbose
      ),
      error = function (e) pStop(who="enkf",conditionMessage(e))
    )

  }
)

##' @rdname kalman
##' @export
setMethod(
  "enkf",
  signature=signature(data="pomp"),
  function (data, Np,
    ..., verbose = getOption("verbose", FALSE)) {

    tryCatch(
      enkf_internal(
        data,
        Np=Np,
        ...,
        verbose=verbose
      ),
      error = function (e) pStop(who="enkf",conditionMessage(e))
    )

  }
)

##' @rdname kalman
##' @export
setMethod(
  "enkf",
  signature=signature(data="kalmand_pomp"),
  function (data, Np,
    ..., verbose = getOption("verbose", FALSE)) {
    if (missing(Np)) Np <- data@Np
    tryCatch(
      enkf_internal(
        as(data,"pomp"),
        Np=Np,
        ...,
        verbose=verbose
      ),
      error = function (e) pStop(who="enkf",conditionMessage(e))
    )

  }
)

## ENSEMBLE KALMAN FILTER (ENKF)

## Ensemble: $X_t\in \mathbb{R}^{m\times q}$
## Prediction mean: $M_t=\langle X \rangle$
## Prediction variance: $V_t=\langle\langle X \rangle\rangle$
## Forecast: $Y_t=emeasure(X_t)$
## Forecast mean: $N_t=\langle Y \rangle$.
## Forecast variance: $S_t=\langle\langle Y \rangle\rangle$
## State/forecast covariance: $W_t=\langle\langle X,Y\rangle\rangle$
## Kalman gain: $K_t = W_t\,S_t^{-1}$
## New observation: $y_t\in \mathbb{R}^{n\times 1}$
## Updated ensemble: $X^u_{t}=X_t + K_t\,(O_t - Y_t)$
## Filter mean: $m_t=\langle X^u_t \rangle = \frac{1}{q} \sum\limits_{i=1}^q x^{u_i}_t$

enkf_internal <- function (object, Np, ..., verbose) {

  verbose <- as.logical(verbose)

  object <- pomp(object,...,verbose=verbose)

  if (undefined(object@rprocess))
    pStop_(paste(sQuote(c("rprocess")),collapse=", ")," is a needed basic component.")

  if (undefined(object@emeasure))
    pStop_(paste(sQuote(c("emeasure")),collapse=", ")," is a needed basic component.")

  if (undefined(object@vmeasure))
    pStop_(paste(sQuote(c("vmeasure")),collapse=", ")," is a needed basic component.")

  Np <- as.integer(Np)
  if (length(Np)>1 || !is.finite(Np) || isTRUE(Np<=0))
    pStop_(sQuote("Np")," should be a single positive integer.")
  
  params <- coef(object)

  t <- time(object)
  tt <- time(object,t0=TRUE)
  ntimes <- length(t)

  y <- obs(object)
  nobs <- nrow(y)

  pompLoad(object,verbose=verbose)
  on.exit(pompUnload(object,verbose=verbose))

  X <- rinit(object,params=params,nsim=Np)
  nvar <- nrow(X)

  filterMeans <- array(dim=c(nvar,ntimes),dimnames=list(name=rownames(X),time=t))
  predMeans <- array(dim=c(nvar,ntimes),dimnames=list(name=rownames(X),time=t))
  forecast <- array(dim=c(nobs,ntimes),dimnames=dimnames(y))
  condlogLik <- numeric(ntimes)

  first <- TRUE

  for (k in seq_len(ntimes)) {
    ## advance ensemble according to state process
    X <- rprocess(object,x0=X,t0=tt[k],times=tt[k+1],params=params,.gnsi=first)
    rn <- rownames(X)
    
    predMeans[,k] <- pm <- rowMeans(X) # prediction mean
    Y <- emeasure(object,x=X,times=tt[k+1],params=params,.gnsi=first)
    ym <- rowMeans(Y)                  # forecast mean
    if (first && nrow(Y) != nobs)
      pStop_("number of observables returned by ",sQuote("emeasure"),
        " does not match data.")

    if (first && (nrow(Y)!=nrow(y) || !all(rownames(Y)==rownames(y)))) {
      y <- y[rownames(Y),]
      rownames(forecast) <- rownames(Y)
    }

    R <- vmeasure(object,x=pm,params=params,times=tt[k+1],.gnsi=first)
    if (first && (nrow(R)!=nrow(y) || !all(rownames(R)==rownames(y)))) {
      pStop_("rownames of matrix returned by ",sQuote("vmeasure"),
        " do not match those returned by ",sQuote("emeasure"),".")
    }

    first <- FALSE

    dn <- dim(R)[c(1L,2L)]
    dim(R) <- dn
    
    sqrtR <- tryCatch(
      t(chol(R)), ## t(sqrtR)%*%sqrtR == R
      error = function (e) {
        pStop_("degenerate ",sQuote("vmeasure"),": ",conditionMessage(e))
      }
    )

    X <- X-pm
    Y <- Y-ym

    dim(X) <- dim(X)[-3L]
    rownames(X) <- rn
    dim(Y) <- dim(Y)[-3L]

    fv <- tcrossprod(Y)/(Np-1L)+R        # forecast variance
    vyx <- tcrossprod(Y,X)/(Np-1L) # forecast/state covariance

    svdS <- svd(fv,nv=0)            # singular value decomposition
    Kt <- svdS$u%*%(crossprod(svdS$u,vyx)/svdS$d) # transpose of Kalman gain
    Ek <- sqrtR%*%matrix(rnorm(n=nobs*Np),nobs,Np) # artificial noise
    resid <- y[,k]-ym

    X <- X+pm+crossprod(Kt,resid+Ek-Y)

    condlogLik[k] <- sum(dnorm(x=crossprod(svdS$u,resid),mean=0,sd=sqrt(svdS$d),log=TRUE))
    filterMeans[,k] <- rowMeans(X)  # filter mean
    forecast[,k] <- ym
  }

  new(
    "kalmand_pomp",
    object,
    Np=Np,
    filter.mean=filterMeans,
    pred.mean=predMeans,
    forecast=forecast,
    cond.logLik=condlogLik,
    loglik=sum(condlogLik)
  )

}

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

##' @rdname kalman
##' @aliases eakf
##' @export
setMethod(
  "eakf",
  signature=signature(data="data.frame"),
  function (data,
    Np, params, rinit, rprocess, emeasure, vmeasure,
    ...,verbose = getOption("verbose", FALSE)) {

    tryCatch(
      eakf_internal(
        data,
        Np=Np,
        params=params,
        rinit=rinit,
        rprocess=rprocess,
        emeasure=emeasure,
        vmeasure=vmeasure,
        ...,
        verbose=verbose
      ),
      error = function (e) pStop(who="eakf",conditionMessage(e))
    )

  }
)

##' @rdname kalman
##' @export
setMethod(
  "eakf",
  signature=signature(data="pomp"),
  function (data,
    Np, ..., verbose = getOption("verbose", FALSE)) {

    tryCatch(
      eakf_internal(
        data,
        Np=Np,
        ...,
        verbose=verbose
      ),
      error = function (e) pStop(who="eakf",conditionMessage(e))
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

eakf_internal <- function (object, Np, ..., verbose) {

  verbose <- as.logical(verbose)

  object <- pomp(object,...,verbose=verbose)

  if (undefined(object@rprocess))
    pStop_(paste(sQuote(c("rprocess")),collapse=", ")," is a needed basic component.")

  if (undefined(object@emeasure))
    pStop_(paste(sQuote(c("emeasure")),collapse=", ")," is a needed basic component.")

  if (undefined(object@vmeasure))
    pStop_(paste(sQuote(c("vmeasure")),collapse=", ")," is a needed basic component.")

  Np <- as.integer(Np)
  if (length(Np)>1 || !is.finite(Np) || isTRUE(Np<=0))
    pStop_(sQuote("Np")," should be a single positive integer.")
  
  params <- coef(object)

  t <- time(object)
  tt <- time(object,t0=TRUE)
  ntimes <- length(t)

  y <- obs(object)
  nobs <- nrow(y)

  pompLoad(object,verbose=verbose)
  on.exit(pompUnload(object,verbose=verbose))

  X <- rinit(object,params=params,nsim=Np)
  nvar <- nrow(X)

  filterMeans <- array(dim=c(nvar,ntimes),
    dimnames=list(name=rownames(X),time=t))
  predMeans <- array(dim=c(nvar,ntimes),
    dimnames=list(name=rownames(X),time=t))
  forecast <- array(dim=c(nobs,ntimes),dimnames=dimnames(y))
  condlogLik <- numeric(ntimes)

  first <- TRUE

  for (k in seq_len(ntimes)) {

    ## advance ensemble according to state process
    X <- rprocess(object,x0=X,t0=tt[k],times=tt[k+1],params=params,.gnsi=first)
    Y <- emeasure(object,x=X,times=tt[k+1],params=params,.gnsi=first)
    if (first && nrow(Y) != nobs)
      pStop_("number of observables returned by ",sQuote("emeasure"),
        " does not match data.")
    if (first && (nrow(Y)!=nrow(y) || !all(rownames(Y)==rownames(y)))) {
      y <- y[rownames(Y),]
      rownames(forecast) <- rownames(Y)
    }
    dn <- dim(Y)[c(1L,2L)]
    dim(Y) <- dn

    dn <- dim(X)[c(1L,2L)]
    nx <- rownames(X)
    dim(X) <- dn
    rownames(X) <- nx

    predMeans[,k] <- pm <- rowMeans(X) # prediction mean
    X <- X-pm

    pv <- tcrossprod(X)/(Np-1L) # prediction variance
    svdV <- svd(pv,nv=0)        # SVD 1: pv = U%*%D%*%t(U)

    Ct <- solve(pv,tcrossprod(X,Y))/(Np-1L) # t(C)

    forecast[,k] <- emeasure(object,x=pm,times=tt[k+1],params=params,.gnsi=FALSE) # forecast (observables)
    resid <- y[,k]-forecast[,k] # forecast error

    R <- vmeasure(object,x=pm,params=params,times=tt[k+1],.gnsi=first)
    if (first && (nrow(R)!=nrow(y) || !all(rownames(R)==rownames(y)))) {
      pStop_("rownames of matrix returned by ",sQuote("vmeasure"),
        " do not match those returned by ",sQuote("emeasure"),".")
    }
    dn <- dim(R)[c(1L,2L)]
    dim(R) <- dn
    first <- FALSE

    Ri <- tryCatch(
      solve(R),
      error = function (e) {
        pStop_("degenerate ",sQuote("vmeasure"),": ",conditionMessage(e))
      }
    )

    ww <- crossprod(svdV$u,Ct)
    w <- crossprod(ww,svdV$d*ww)+R # forecast variance
    svdW <- svd(w,nv=0) # SVD 2:  needed for log likelihood

    condlogLik[k] <- sum(dnorm(x=crossprod(svdW$u,resid),mean=0,
      sd=sqrt(svdW$d),log=TRUE))

    u <- sqrt(svdV$d)*ww
    u <- tcrossprod(u%*%Ri,u)
    svdU <- svd(u,nv=0) # SVD 3

    ## adjustment
    b <- svdV$u%*%(sqrt(svdV$d)*svdU$u)%*%
      (1/sqrt(1+svdU$d)/sqrt(svdV$d)*t(svdV$u))

    K <- tcrossprod(b%*%pv,b)%*%Ct%*%Ri # Kalman gain

    filterMeans[,k] <- fm <- pm+K%*%resid # filter mean

    X[,] <- b%*%X+fm[,]
  }

  new(
    "kalmand_pomp",
    object,
    Np=Np,
    filter.mean=filterMeans,
    pred.mean=predMeans,
    forecast=forecast,
    cond.logLik=condlogLik,
    loglik=sum(condlogLik)
  )
}
