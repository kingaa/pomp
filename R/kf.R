##' Kalman filter
##'
##' The basic Kalman filter for multivariate, linear, Gaussian processes.
##'
##' If the latent state is \eqn{X}, the observed variable is \eqn{Y},
##' \eqn{X_t \in R^m}, \eqn{Y_t \in R^n}, and
##' \deqn{X_t \sim \mathrm{MVN}(A X_{t-1}, Q)}{X_t ~ MVN(A X_{t-1}, Q)}
##' \deqn{Y_t \sim \mathrm{MVN}(C X_t, R)}{Y_t ~ MVN(C X_t, R)}
##' where \eqn{\mathrm{MVN}(M,V)} denotes the multivariate normal distribution with mean \eqn{M} and variance \eqn{V}.
##' Then the Kalman filter computes the exact likelihood of \eqn{Y}
##' given \eqn{A}, \eqn{C}, \eqn{Q}, and \eqn{R}.
##'
##' @rdname kf
##' @concept Kalman filter
##' @seealso \code{\link{enkf}}, \code{\link{eakf}}
##'
##' @param object a pomp object containing data;
##' @param X0 length-m vector containing initial state.
##' This is assumed known without uncertainty.
##' @param A \eqn{m\times m}{m x m} latent state-process transition matrix.
##' \eqn{E[X_{t+1} | X_t] = A.X_t}.
##' @param Q \eqn{m\times m}{m x m} latent state-process covariance matrix.
##' \eqn{Var[X_{t+1} | X_t] = Q}
##' @param C \eqn{n\times m}{n x m} link matrix.
##' \eqn{E[Y_t | X_t] = C.X_t}.
##' @param R \eqn{n\times n}{n x n} observation process covariance matrix.
##' \eqn{Var[Y_t | X_t] = R}
##' @param tol numeric;
##' the tolerance to be used in computing matrix pseudoinverses via singular-value decomposition.
##' Singular values smaller than \code{tol} are set to zero.
##'
##' @example examples/kf.R
##'
##' @return
##' A named list containing the following elements:
##' \describe{
##'   \item{object}{the \sQuote{pomp} object}
##'   \item{A, Q, C, R}{as in the call}
##'   \item{filter.mean}{\eqn{E[X_t|y^*_1,\dots,y^*_t]}}
##'   \item{pred.mean}{\eqn{E[X_t|y^*_1,\dots,y^*_{t-1}]}}
##'   \item{forecast}{\eqn{E[Y_t|y^*_1,\dots,y^*_{t-1}]}}
##'   \item{cond.logLik}{\eqn{f(y^*_t|y^*_1,\dots,y^*_{t-1})}}
##'   \item{logLik}{\eqn{f(y^*_1,\dots,y^*_T)}}
##' }
##
## X(t) ~ MultivariateNormal(A X(t-1), Q)
## Y(t) ~ MultivariateNormal(C X(t), R)
##
## FM,FV = filter mean, variance
## PM,PV = prediction mean, variance
## YM,YV = forecast mean, variance
## r = residual
## K = Kalman gain
##
## PM(t) = A FM(t-1)
## PV(t) = A FV(t-1) A^T + Q
## YM(t) = C PM(t)
## r(t) = y*(t) - YM(t)
## YV(t) = C PV(t) C^T + R
## FV(t)^(-1) = PV(t)^(-1) + C^T R^(-1) C
## K(t) = FV(t) C^T R^(-1)
## FM(t) = PM(t) + K(t) r(t)
##
##' @export
kalmanFilter <- function (
  object, X0 = rinit(object),
  A, Q, C, R, tol = 1e-6
) {

  t <- time(object,t0=FALSE)
  y <- obs(object)
  N <- length(t)     # number of observations
  nvar <- length(X0) # number of state variables
  nobs <- nrow(y)    # number of observables

  filterMeans <- array(
    dim=c(nvar,N),
    dimnames=list(name=names(X0),time=NULL)
  )
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

    ## first SVD
    s1 <- svd_pseudoinv(pv,tol=tol)
    pvi <- s1$u%*%(t(s1$u)/s1$d)       # prediction precision
    fvi <- pvi+cric                    # filter precision

    forecast[,k] <- ym <- C%*%pm
    resid <- y[,k]-ym                  # forecast error
    yv <- tcrossprod(C%*%pv,C)+R       # forecast variance

    ## second SVD
    s2 <- svd_pseudoinv(yv,tol=tol)

    condlogLik[k] <- sum(
      dnorm(x=crossprod(s2$u,resid),mean=0,sd=sqrt(s2$d),log=TRUE)
    )

    ## third SVD
    s3 <- svd_pseudoinv(fvi,tol=tol)
    fv <- s3$u%*%(t(s3$u)/s3$d)        # filter variance
    K <- fv%*%crossprod(C,ri)          # Kalman gain

    filterMeans[,k] <- fm <- pm+K%*%resid  # filter mean
  }

  list(
    object=object,
    A=A,
    Q=Q,
    C=C,
    R=R,
    filter.mean=filterMeans,
    pred.mean=predMeans,
    forecast=forecast,
    cond.logLik=condlogLik,
    logLik=sum(condlogLik)
  )

}

svd_pseudoinv <- function (A, tol) {
  s <- svd(A,nv=0,nu=0)
  v <- s$d/max(s$d) > tol
  s <- svd(A,nv=0,nu=sum(v))
  s$d <- s$d[v]
  s
}
