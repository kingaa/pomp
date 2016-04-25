setClass(
    "kalmand.pomp",
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

setMethod(
    "enkf",
    signature=signature(object="pomp"),
    function (object, params, Np, h, R,
              verbose = getOption("verbose"),
              ...) {
        if (missing(params)) params <- coef(object)
        enkf.internal(
            object=object,
            params=params,
            h=h,
            R=R,
            Np=Np,
            verbose=verbose,
            ...
        )
    }
)

enkf.internal <- function (object, params, h, R, Np, verbose) {
    Np <- as.integer(Np)
    R <- as.matrix(R)
    verbose <- as.logical(verbose)
    
    t <- time(object)
    tt <- time(object,t0=TRUE)
    ntimes <- length(t)

    y <- obs(object)
    nobs <- nrow(y)

    Y <- array(dim=c(nobs,Np),dimnames=list(variable=rownames(y),rep=NULL))
    X <- init.state(object,params=params,nsim=Np)
    nvar <- nrow(X)

    filterMeans <- array(dim=c(nvar,ntimes),dimnames=list(variable=rownames(X),time=t))
    predMeans <- array(dim=c(nvar,ntimes),dimnames=list(variable=rownames(X),time=t))
    forecast <- array(dim=c(nobs,ntimes),dimnames=dimnames(y))
    condlogLik <- numeric(ntimes)

    sqrtR <- try(t(chol(R)))            # t(sqrtR)%*%sqrtR == R
    if (inherits(sqrtR,"try-error"))
        stop("error in ",sQuote("enkf"),": degenerate ",sQuote("R"),call.=FALSE)

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

    new("kalmand.pomp",object,Np=Np,
        filter.mean=filterMeans,
        pred.mean=predMeans,
        forecast=forecast,
        cond.loglik=condlogLik,
        loglik=sum(condlogLik))
}

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

setMethod(
    "eakf",
    signature=signature(object="pomp"),
    function (object, params, Np, C, R,
              verbose = getOption("verbose"),
              ...) {
        if (missing(params)) params <- coef(object)
        eakf.internal(
            object=object,
            params=params,
            C=C,
            R=R,
            Np=Np,
            verbose=verbose,
            ...
        )
    }
)

eakf.internal <- function (object, params, C, R, Np, verbose) {
    Np <- as.integer(Np)
    R <- as.matrix(R)
    verbose <- as.logical(verbose)

    t <- time(object)
    tt <- time(object,t0=TRUE)

    y <- obs(object)

    X <- init.state(object,params=params,nsim=Np)

    nvar <- nrow(X)
    ntimes <- length(t)
    nobs <- nrow(y)

    filterMeans <- array(dim=c(nvar,ntimes),dimnames=list(variable=rownames(X),time=t))
    predMeans <- array(dim=c(nvar,ntimes),dimnames=list(variable=rownames(X),time=t))
    forecast <- array(dim=c(nobs,ntimes),dimnames=dimnames(y))
    condlogLik <- numeric(ntimes)

    ri <- try(solve(R))
    if (inherits(ri,"try-error"))
        stop("error in ",sQuote("eakf"),": degenerate ",sQuote("R"),call.=FALSE)

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

        condlogLik[k] <- sum(dnorm(x=crossprod(svdW$u,resid),mean=0,sd=sqrt(svdW$d),log=TRUE))

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

    new("kalmand.pomp",object,Np=Np,
        filter.mean=filterMeans,
        pred.mean=predMeans,
        forecast=forecast,
        cond.loglik=condlogLik,
        loglik=sum(condlogLik))
}

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
    condlogLik[k] <- sum(dnorm(x=crossprod(svdW$u,resid),mean=0,sd=sqrt(svdW$d),log=TRUE))
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
