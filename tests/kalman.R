options(digits=3)
png(filename="kalman-%02d.png",res=100)

library(pomp)
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(mvtnorm)
})
set.seed(1638125322)

t <- seq(1,20)

C <- matrix(c(
  1, 0,  0, 0,
  0, 1, -1, 0),
  2,4,byrow=TRUE)

dimX <- ncol(C)
dimY <- nrow(C)

dimnames(C) <- list(
  paste0("y",seq_len(dimY)),
  paste0("x",seq_len(dimX)))

R <- matrix(c(1,0.3,0.3,1),nrow=dimY,dimnames=list(rownames(C),rownames(C)))

A <- matrix(c(
  -0.677822, -0.169411,  0.420662,  0.523571,
  2.87451,  -0.323604, -0.489533, -0.806087,
  -1.36617,  -0.592326,  0.567114,  0.345142,
  -0.807978, -0.163305,  0.668037,  0.468286),
  nrow=dimX,ncol=dimX,byrow=TRUE,
  dimnames=list(colnames(C),colnames(C)))

Q <- crossprod(matrix(rnorm(n=dimX*dimX,sd=1/4),4,4))
dimnames(Q) <- dimnames(A)

X0 <- setNames(rnorm(dimX),colnames(C))

N <- length(t)
x <- array(dim=c(dimX,N),dimnames=list(variable=colnames(C),time=t))
y <- array(dim=c(dimY,N),dimnames=list(variable=rownames(C),time=t))

xx <- X0
sqrtQ <- t(chol(Q))
sqrtR <- t(chol(R))
for (k in seq_along(t)) {
  x[,k] <- xx <- A %*% xx + sqrtQ %*% rnorm(n=dimX)
  y[,k] <- C %*% xx + sqrtR %*% rnorm(n=dimY)
}

y |>
  t() |>
  as.data.frame() |>
  pfilter(
    times=t,t0=0,
    A=A,C=C,R=R,sqrtQ=sqrtQ,sqrtR=sqrtR,X0=X0,
    rprocess=discrete_time(
      step.fun=function(x1,x2,x3,x4,delta.t,...){
        x <- c(x1,x2,x3,x4)
        x <- A%*%x+sqrtQ%*%rnorm(n=ncol(A))
        setNames(x,c("x1","x2","x3","x4"))
      },
      delta.t=1),
    emeasure=function(x1,x2,x3,x4,C,...){
      ex <- C%*%c(x1,x2,x3,x4)
      dim(ex) <- NULL
      setNames(ex,rownames(C))
    },
    vmeasure=function(R, ...){
      R
    },
    rmeasure=function(x1,x2,x3,x4,C,sqrtR,...){
      x <- c(x1,x2,x3,x4)
      C%*%x+sqrtR%*%rnorm(n=nrow(C))
    },
    dmeasure=function(x1,x2,x3,x4,y1,y2,log,C,R,...){
      x <- c(x1,x2,x3,x4)
      y <- c(y1,y2)
      dmvnorm(x=t(y-C%*%x),sigma=R,log=log)
    },
    rinit=function(params,t0,X0,...){
      X0
    },
    params=c(),
    Np=1000,
    filter.mean=TRUE
  ) -> pf

freeze(
  seed=1286628545,
  {
    enkf <- enkf(pf,Np=1000)
    eakf <- eakf(pf,Np=1000)
    kf <- kalmanFilter(pf,A=A,Q=Q,C=C,R=R,tol=1e-9)
    stopifnot(
      all.equal(
        c(kf$logLik,logLik(pf),logLik(enkf),logLik(eakf)),
        c(-67.0,-67.1,-66.7,-67.1),
        tolerance=1e-3
      )
    )
  }
)

try(enkf(pf))
try(enkf(pf,Np=c(100,200)))
try(enkf(pf,Np=-10))
try(enkf(pf,Np="10b"))
try(enkf(pf,Np=100,emeasure=NULL))
try(enkf(enkf,Np=100,vmeasure=NULL))
try(enkf(enkf,C=matrix(1,2,4,dimnames=list(c("y2","y1"),NULL)),Np=100))
try(enkf(enkf,C=matrix(1,3,4,dimnames=list(c("a","b","c"),NULL)),Np=100))
try(enkf(enkf,R=matrix(1,1,1,dimnames=list(c("y1"),NULL)),Np=100))
try(enkf(enkf,R=matrix(1,2,2,dimnames=list(c("y2","y1"),NULL)),Np=100))
try(enkf(enkf,R=matrix(1,2,2,dimnames=list(c("a","b"),NULL)),Np=100))
try(enkf(enkf,R=matrix(1,3,3,dimnames=list(c("a","b","c"),NULL)),Np=100))
enkf(enkf)
enkf(enkf,Np=100)

try(eakf(enkf,Np=c(100,200)))
try(eakf(enkf,Np=-10))
try(eakf(enkf,Np="10b"))
try(eakf(enkf,Np=100,emeasure=NULL))
try(eakf(enkf,Np=100,vmeasure=NULL))
try(eakf(enkf,C=matrix(1,2,4,dimnames=list(c("y2","y1"),NULL)),Np=100))
try(eakf(enkf,C=matrix(1,3,4,dimnames=list(c("a","b","c"),NULL)),Np=100))
try(eakf(enkf,R=matrix(1,1,1,dimnames=list(c("y1"),NULL)),Np=100))
try(eakf(enkf,R=matrix(1,2,2,dimnames=list(c("y2","y1"),NULL)),Np=100))
try(eakf(enkf,R=matrix(1,2,2,dimnames=list(c("a","b"),NULL)),Np=100))
try(eakf(enkf,R=matrix(1,3,3,dimnames=list(c("a","b","c"),NULL)),Np=100))

invisible(enkf(pf,Np=1000,params=as.list(coef(pf))))
invisible(eakf(pf,Np=1000,params=as.list(coef(pf))))

enkf |>
  as.data.frame() |>
  pivot_longer(cols=-time) |>
  group_by(name) |>
  summarize(n=length(value))
eakf |>
  as.data.frame() |>
  pivot_longer(cols=-time) |>
  group_by(name) |>
  summarize(n=length(value))

enkf |>
  forecast(format="d") |>
  ggplot(aes(x=time,y=value,group=variable,color=variable))+
  geom_line()+theme_bw()+
  labs(title="EnKF forecasts")

eakf |>
  forecast(vars=c("y1","y2"),format="d") |>
  ggplot(aes(x=time,y=value,group=variable,color=variable))+
  geom_line()+theme_bw()+
  labs(title="EAKF forecasts")

try({
  R <- matrix(c(1,0,1,0),2,2)
  rownames(R) <- rownames(C)
  enkf(pf,Np=1000,R=R)
})

try({
  eakf(pf,Np=1000,R=R)
})

try(enkf())
try(enkf("bob"))
try(eakf())
try(eakf("bob"))

dev.off()
