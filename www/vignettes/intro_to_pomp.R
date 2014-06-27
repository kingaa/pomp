
## ----include=FALSE-------------------------------------------------------
opts_chunk$set(
               echo=TRUE,results='markup',
               progress=TRUE,prompt=FALSE,tidy=FALSE,highlight=FALSE,
               print=FALSE,keep.source=TRUE,comment='##',
               size='normalsize',background="#FFFFFF",
               warning=FALSE,message=FALSE,error=FALSE,
               dev='png',
               fig.path='figure/intro-',fig.lp="fig:",
               fig.align='left',fig.show='asis',
               fig.height=6,fig.width=8,
               dpi=150,
               dev.args=list(
                 bg='transparent',
                 pointsize=10
                 )
               )



## ----set-opts,echo=F,results='hide'--------------------------------------
set.seed(5384959L)


## ----gompertz-proc-sim-def-----------------------------------------------
require(pomp)

gompertz.proc.sim <- function (x, t, params, delta.t, ...) {
  ## unpack the parameters:
  r <- params["r"]
  K <- params["K"]
  sigma <- params["sigma"]
  ## the state at time t:
  X <- x["X"]
  ## generate a log-normal random variable:
  eps <- exp(rnorm(n=1,mean=0,sd=sigma))
  ## compute the state at time t+delta.t:
  S <- exp(-r*delta.t)
  xnew <- c(X=unname(K^(1-S)*X^S*eps))
  return(xnew)
}



## ----gompertz-meas-sim-def-----------------------------------------------
gompertz.meas.sim <- function (x, t, params, ...) {
  ## unpack the parameters:
  tau <- params["tau"]
  ## state at time t:
  X <- x["X"]
  ## generate a simulated observation:
  y <- c(Y=unname(rlnorm(n=1,meanlog=log(X),sd=tau)))
  return(y)
}



## ----first-pomp-construction,eval=F--------------------------------------
## gompertz <- pomp(
##                  data=data.frame(
##                    time=1:100,
##                    Y=NA
##                    ),
##                  times="time",
##                  rprocess=discrete.time.sim(
##                    step.fun=gompertz.proc.sim,
##                    delta.t=1
##                    ),
##                  rmeasure=gompertz.meas.sim,
##                  t0=0
##                  )
## 


## ----set-params----------------------------------------------------------
theta <- c(r=0.1,K=1,sigma=0.1,tau=0.1,X.0=1)


## ----gompertz-first-simulation,eval=F------------------------------------
## gompertz <- simulate(gompertz,params=theta)

## ----gompertz-get-data,eval=T,echo=F,results='hide'----------------------
pompExample(gompertz)
dat <- as.data.frame(gompertz)
gompertz <- pomp(
                 data=dat[c("time","Y")],
                 times="time",
                 rprocess=discrete.time.sim(
                   step.fun=gompertz.proc.sim,
                   delta.t=1
                   ),
                 rmeasure=gompertz.meas.sim,
                 t0=0
                 )
coef(gompertz) <- theta


## ----gompertz-plot,echo=F------------------------------------------------
plot(gompertz,variables=c("Y"))


## ----second-pomp-construction--------------------------------------------
gompertz.meas.dens <- function (y, x, t, params, log, ...) {
  ## unpack the parameters:
  tau <- params["tau"]
  ## state at time t:
  X <- x["X"]
  ## observation at time t:
  Y <- y["Y"]
  ## compute the likelihood of Y|X,tau
  f <- dlnorm(x=Y,meanlog=log(X),sdlog=tau,log=log)
  return(f)
}

gompertz <- pomp(
                 gompertz,
                 dmeasure=gompertz.meas.dens
                 )


## ----gompertz-pfilter-truth,eval=F---------------------------------------
## pf <- pfilter(gompertz,params=theta,Np=1000)
## loglik.truth <- logLik(pf)
## loglik.truth

## ----gompertz-pfilter-truth-eval,echo=F----------------------------------
set.seed(457645443L)
pf <- pfilter(gompertz,params=theta,Np=1000)
loglik.truth <- logLik(pf)
loglik.truth


## ----gompertz-pfilter-truth-alt1,eval=F----------------------------------
## pf <- pfilter(gompertz,params=coef(gompertz),Np=1000)


## ----gompertz-pfilter-truth-alt2,eval=F----------------------------------
## pf <- pfilter(gompertz,Np=1000)


## ----gompertz-pfilter-guess,eval=F---------------------------------------
## theta.true <- coef(gompertz)
## theta.guess <- theta.true
## theta.guess[c("r","K","sigma")] <- 1.5*theta.true[c("r","K","sigma")]
## pf <- pfilter(gompertz,params=theta.guess,Np=1000)
## loglik.guess <- logLik(pf)

## ----gompertz-pfilter-guess-eval,echo=F----------------------------------
binary.file <- "gompertz-pfilter-guess.rda"
if (file.exists(binary.file)) {
  load(binary.file)
} else {
set.seed(457645443L)
theta.true <- coef(gompertz)
theta.guess <- theta.true
theta.guess[c("r","K","sigma")] <- 1.5*theta.true[c("r","K","sigma")]
pf <- pfilter(gompertz,params=theta.guess,Np=1000)
loglik.guess <- logLik(pf)
save(theta.true,theta.guess,loglik.guess,file=binary.file,compress='xz')
}  


## ----kalman-filter-def---------------------------------------------------
kalman.filter <- function (Y, X0, r, K, sigma, tau) {
  ntimes <- length(Y)
  sigma.sq <- sigma^2
  tau.sq <- tau^2
  cond.loglik <- numeric(ntimes)
  filter.mean <- numeric(ntimes)
  pred.mean <- numeric(ntimes)
  pred.var <- numeric(ntimes)
  m <- log(X0)
  v <- 0
  S <- exp(-r)
  for (k in seq_len(ntimes)) {
    pred.mean[k] <- M <- (1-S)*log(K) + S*m
    pred.var[k] <- V <- S*v*S+sigma.sq
    q <- V+tau.sq
    r <- log(Y[k])-M
    cond.loglik[k] <- dnorm(x=log(Y[k]),mean=M,sd=sqrt(q),log=TRUE)-log(Y[k])
    q <- 1/V+1/tau.sq
    filter.mean[k] <- m <- (log(Y[k])/tau.sq+M/V)/q
    v <- 1/q
  }
  list(
       pred.mean=pred.mean,
       pred.var=pred.var,
       filter.mean=filter.mean,
       cond.loglik=cond.loglik,
       loglik=sum(cond.loglik)
       )
}


## ----kalman-filter-run---------------------------------------------------
y <- obs(gompertz)
x0 <- init.state(gompertz)
r <- coef(gompertz,"r")
K <- coef(gompertz,"K")
sigma <- coef(gompertz,"sigma")
tau <- coef(gompertz,"tau")
kf <- kalman.filter(y,x0,r,K,sigma,tau)

## ----kalman-likelihood-correction,echo=F---------------------------------
loglik.kalman <- kf$loglik


## ----eval=F--------------------------------------------------------------
## as(gompertz,"data.frame")


## ----eval=F--------------------------------------------------------------
## obs(gompertz)
## obs(gompertz,"Y")
## time(gompertz)


## ----eval=F--------------------------------------------------------------
## time(gompertz) <- 1:10


## ----eval=F--------------------------------------------------------------
## timezero(gompertz)
## timezero(gompertz) <- -10


## ----eval=F--------------------------------------------------------------
## time(gompertz,t0=TRUE)
## time(gompertz,t0=T) <- seq(from=0,to=10,by=1)


## ----eval=F--------------------------------------------------------------
## window(gompertz,start=3,end=20)


## ----eval=F--------------------------------------------------------------
## coef(gompertz)
## coef(gompertz,c("sigma","tau")) <- c(1,0)


## ----eval=F--------------------------------------------------------------
## states(gompertz)
## states(gompertz,"X")


## ----snippet-gomp-pomp,results='hide'------------------------------------
gomp2 <- pomp(
              data=subset(as(gompertz,"data.frame"),select=c(time,Y)),
              times="time",
              t0=0,
              rmeasure=Csnippet('
   Y = rlnorm(log(X),tau);
'),
              dmeasure=Csnippet('
   lik = dlnorm(Y,log(X),tau,give_log);
'),
              rprocess=discrete.time.sim(
                step.fun=Csnippet('
  double S = exp(-r*dt);
  double eps = rlnorm(0,sigma);
  X = pow(K,(1-S))*pow(X,S)*eps;
'),
                delta.t=1
                ),
              paramnames=c("sigma","tau","r","K"),
              statenames="X",
              params=theta
              )



## ----gompertz-perform,eval=F,echo=T--------------------------------------
## 
## tic <- Sys.time()
## sim1 <- simulate(gompertz,nsim=1000,seed=5676868L,obs=TRUE)
## toc <- Sys.time()
## g1sim <- toc-tic
## 
## tic <- Sys.time()
## sim2 <- simulate(gomp2,nsim=1000,seed=5676868L,obs=TRUE)
## toc <- Sys.time()
## g2sim <- toc-tic
## 
## stopifnot(all.equal(sim1,sim2))
## 
## tic <- Sys.time()
## pf1 <- pfilter(gompertz,Np=10000,seed=5676868L)
## toc <- Sys.time()
## g1pf <- toc-tic
## 
## tic <- Sys.time()
## pf2 <- pfilter(gomp2,Np=10000,seed=5676868L)
## toc <- Sys.time()
## g2pf <- toc-tic
## 
## stopifnot(all.equal(logLik(pf1),logLik(pf2)))
## 

## ----gompertz-perform-eval,eval=T,echo=F---------------------------------
binary.file <- "gompertz-performance.rda"
if (file.exists(binary.file)) {
  load(binary.file)
} else {
set.seed(457645443L)

tic <- Sys.time()
sim1 <- simulate(gompertz,nsim=1000,seed=5676868L,obs=TRUE)
toc <- Sys.time()
g1sim <- toc-tic

tic <- Sys.time()
sim2 <- simulate(gomp2,nsim=1000,seed=5676868L,obs=TRUE)
toc <- Sys.time()
g2sim <- toc-tic

stopifnot(all.equal(sim1,sim2))

tic <- Sys.time()
pf1 <- pfilter(gompertz,Np=10000,seed=5676868L)
toc <- Sys.time()
g1pf <- toc-tic

tic <- Sys.time()
pf2 <- pfilter(gomp2,Np=10000,seed=5676868L)
toc <- Sys.time()
g2pf <- toc-tic

stopifnot(all.equal(logLik(pf1),logLik(pf2)))

save(g1sim,g2sim,g1pf,g2pf,file=binary.file,compress='xz')
}  


## ----eval=F--------------------------------------------------------------
## demo(gompertz)


## ----eval=F--------------------------------------------------------------
## system.file("examples",package="pomp")


## ----loggomp-pomp-construction,eval=T------------------------------------
gompertz <- pomp(
                 gompertz,
                 parameter.transform=function(params,...){
                   exp(params)
                 },
                 parameter.inv.transform=function(params,...){
                   log(params)
                 }
                 )


## ----eval=T--------------------------------------------------------------
coef(gompertz) <- c(r=0.1,K=1,tau=0.1,sigma=0.1,X.0=1)


## ----eval=T--------------------------------------------------------------
coef(gompertz)


## ----eval=T--------------------------------------------------------------
coef(gompertz,transform=TRUE) <- c(r=log(0.1),K=0,tau=log(0.1),
                sigma=log(0.1),X.0=0)


## ----eval=T--------------------------------------------------------------
coef(gompertz,transform=TRUE)


## ----eval=T--------------------------------------------------------------
coef(gompertz)


## ----par-trans-inverse-test,results='markup'-----------------------------
# use parameter.inv.transform:
theta <- coef(gompertz,transform=TRUE)  
## theta is on the estimation scale
g2 <- gompertz
## use parameter.transform:
coef(g2,transform=TRUE) <- theta
## compare the internal-scale representations:
all.equal(coef(gompertz),coef(g2))



## ----echo=F,results='hide'-----------------------------------------------
pompExample(gompertz)
theta <- coef(gompertz)
theta.true <- theta

## ----gompertz-multi-mif-calc,eval=F,echo=T-------------------------------
## estpars <- c("r","sigma","tau")
## mf <- foreach(i=1:10,
##               .inorder=FALSE,
##               .options.multicore=list(set.seed=TRUE)
##               ) %dopar%
## {
##   theta.guess <- theta.true
##   theta.guess[estpars] <- rlnorm(
##                                  n=length(estpars),
##                                  meanlog=log(theta.guess[estpars]),
##                                  sdlog=1
##                                  )
##   m1 <- mif(
##             gompertz,
##             Nmif=50,
##             start=theta.guess,
##             transform=TRUE,
##             rw.sd=c(r=0.02,sigma=0.02,tau=0.05),
##             Np=2000,
##             var.factor=4,
##             ic.lag=10,
##             cooling.type="geometric",
##             cooling.fraction=0.95
##             )
##   m1 <- continue(m1,Nmif=50,cooling.fraction=0.8)
##   m1 <- continue(m1,Nmif=50,cooling.fraction=0.6)
##   m1 <- continue(m1,Nmif=50,cooling.fraction=0.2)
##   ll <- replicate(n=10,logLik(pfilter(m1,Np=10000)))
##   list(mif=m1,ll=ll)
## }


## ----gompertz-post-mif,eval=F,echo=F-------------------------------------
## theta.true <- coef(gompertz)
## loglik.true <- replicate(n=10,logLik(pfilter(gompertz,Np=10000)))
## loglik.true <- logmeanexp(loglik.true,se=TRUE)
## theta.mif <- t(sapply(mf,function(x)coef(x$mif)))
## loglik.mif <- t(sapply(mf,function(x)logmeanexp(x$ll,se=TRUE)))
## best <- which.max(loglik.mif[,1])
## theta.mif <- theta.mif[best,]
## loglik.mif <- loglik.mif[best,]


## ----gompertz-multi-mif-eval,echo=F,results='hide'-----------------------
binary.file <- "gompertz-multi-mif.rda"
if (file.exists(binary.file)) {
  load(binary.file)
} else {

  require(foreach)
  require(doMC)
  registerDoMC()
  
  save.seed <- .Random.seed
  set.seed(334388458L,kind="L'Ecuyer")

  tic <- Sys.time()
estpars <- c("r","sigma","tau")
mf <- foreach(i=1:10,
              .inorder=FALSE,
              .options.multicore=list(set.seed=TRUE)
              ) %dopar%
{
  theta.guess <- theta.true
  theta.guess[estpars] <- rlnorm(
                                 n=length(estpars),
                                 meanlog=log(theta.guess[estpars]),
                                 sdlog=1
                                 )
  m1 <- mif(
            gompertz,
            Nmif=50,
            start=theta.guess,
            transform=TRUE,
            rw.sd=c(r=0.02,sigma=0.02,tau=0.05),
            Np=2000,
            var.factor=4,
            ic.lag=10,
            cooling.type="geometric",
            cooling.fraction=0.95
            )
  m1 <- continue(m1,Nmif=50,cooling.fraction=0.8)
  m1 <- continue(m1,Nmif=50,cooling.fraction=0.6)
  m1 <- continue(m1,Nmif=50,cooling.fraction=0.2)
  ll <- replicate(n=10,logLik(pfilter(m1,Np=10000)))
  list(mif=m1,ll=ll)
}
  toc <- Sys.time()
  etime <- toc-tic
  .Random.seed <<- save.seed

theta.true <- coef(gompertz)
loglik.true <- replicate(n=10,logLik(pfilter(gompertz,Np=10000)))
loglik.true <- logmeanexp(loglik.true,se=TRUE)
theta.mif <- t(sapply(mf,function(x)coef(x$mif)))
loglik.mif <- t(sapply(mf,function(x)logmeanexp(x$ll,se=TRUE)))
best <- which.max(loglik.mif[,1])
theta.mif <- theta.mif[best,]
loglik.mif <- loglik.mif[best,]
  save(
       mf,estpars,
       theta.mif,theta.true,
       loglik.mif,loglik.true,
       etime,
       file=binary.file,
       compress="xz"
       )
}
rbind(
      mle=c(signif(theta.mif[estpars],3),loglik=round(loglik.mif,2)),
      truth=c(signif(theta.true[estpars],3),loglik=round(loglik.true,2))
      ) -> results.table


## ----eval=F--------------------------------------------------------------
## theta.true <- coef(gompertz)
## loglik.true <- replicate(n=10,logLik(pfilter(gompertz,Np=10000)))
## loglik.true <- logmeanexp(loglik.true,se=TRUE)
## theta.mif <- t(sapply(mf,function(x)coef(x$mif)))
## loglik.mif <- t(sapply(mf,function(x)logmeanexp(x$ll,se=TRUE)))
## best <- which.max(loglik.mif[,1])
## theta.mif <- theta.mif[best,]
## loglik.mif <- loglik.mif[best,]


## ----multi-mif-plot,echo=F,eval=F----------------------------------------
## op <- par(mfrow=c(4,1),mar=c(3,3,0,0),mgp=c(2,1,0),bty='l')
## loglik <- sapply(mf,function(x)conv.rec(x$mif,"loglik"))
## log.r <- sapply(mf,function(x)conv.rec(x$mif,"r"))
## log.sigma <- sapply(mf,function(x)conv.rec(x$mif,"sigma"))
## log.tau <- sapply(mf,function(x)conv.rec(x$mif,"tau"))
## matplot(loglik,type='l',lty=1,xlab="",ylab=expression(log~L),xaxt='n',ylim=max(loglik,na.rm=T)+c(-12,3))
## matplot(log.r,type='l',lty=1,xlab="",ylab=expression(log~r),xaxt='n')
## matplot(log.sigma,type='l',lty=1,xlab="",ylab=expression(log~sigma),xaxt='n')
## matplot(log.tau,type='l',lty=1,xlab="MIF iteration",ylab=expression(log~tau))
## par(op)


## ----mif-plot,echo=F-----------------------------------------------------
op <- par(mfrow=c(4,1),mar=c(3,3,0,0),mgp=c(2,1,0),bty='l')
loglik <- sapply(mf,function(x)conv.rec(x$mif,"loglik"))
log.r <- sapply(mf,function(x)conv.rec(x$mif,"r"))
log.sigma <- sapply(mf,function(x)conv.rec(x$mif,"sigma"))
log.tau <- sapply(mf,function(x)conv.rec(x$mif,"tau"))
matplot(loglik,type='l',lty=1,xlab="",ylab=expression(log~L),xaxt='n',ylim=max(loglik,na.rm=T)+c(-12,3))
matplot(log.r,type='l',lty=1,xlab="",ylab=expression(log~r),xaxt='n')
matplot(log.sigma,type='l',lty=1,xlab="",ylab=expression(log~sigma),xaxt='n')
matplot(log.tau,type='l',lty=1,xlab="MIF iteration",ylab=expression(log~tau))
par(op)


## ----first-mif-results-table,echo=F--------------------------------------
print(results.table)


## ----gompertz-skeleton-def,echo=T----------------------------------------
gompertz.skel <- function (x, t, params, ...) {
  r <- params["r"]
  K <- params["K"]
  X <- x["X"]
  S <- exp(-r)
  xnew <- c(X=unname(K^(1-S)*X^S))
  return(xnew)
}


## ----gomp3-pomp----------------------------------------------------------
gomp3 <- simulate(
                  pomp(
                       gompertz,
                       skeleton=gompertz.skel,
                       skeleton.type='map'
                       ),
                  params=c(
                    X.0=0.1,r=0.1,tau=0.05,sigma=0.05,K=1
                    ),
                  seed=88737400L
                  )


## ----gompertz-trajmatch-calc,eval=F--------------------------------------
## tm <- traj.match(
##                  gomp3,
##                  start=coef(gomp3),
##                  transform=TRUE,
##                  est=c("r","K","tau","X.0"),
##                  method="Nelder-Mead",
##                  maxit=1000,
##                  reltol=1e-8
##                  )

## ----gompertz-trajmatch-eval,echo=F,eval=T-------------------------------
binary.file <- "gompertz-trajmatch.rda"
if (file.exists(binary.file)) {
  load(binary.file)
} else {
tm <- traj.match(
                 gomp3,
                 start=coef(gomp3),
                 transform=TRUE,
                 est=c("r","K","tau","X.0"),
                 method="Nelder-Mead",
                 maxit=1000,
                 reltol=1e-8
                 )
  save(tm,file=binary.file,compress="xz")
}


## ----trajmatch-plot,echo=F,eval=T,fig.height=4,fig.width=6---------------
op <- par(mfrow=c(1,1),mar=c(3,3,0,0),mgp=c(2,1,0),bty='l')
plot(time(tm),obs(tm,"Y"),xlab="time",ylab=expression(X,Y),type='o')
lines(time(tm),states(tm,"X"),lwd=2)
par(op)


## ----ricker-map-defn-----------------------------------------------------
ricker.sim <- function (x, t, params, delta.t, ...) {
  e <- rnorm(n=1,mean=0,sd=params["sigma"]) 
  setNames(
           c(
             params["r"]*x["N"]*exp(-x["N"]+e),
             e
             ),
           c("N","e")
           )
}


## ----ricker-sim-C-def----------------------------------------------------
ricker.sim <- Csnippet('
  e = rnorm(0,sigma);
  N = r*N*exp(-N+e);
')


## ----ricker-pomp,results='hide'------------------------------------------
ricker <- pomp(
               data=data.frame(time=seq(0,50,by=1),y=NA),
               times="time",
               t0=0,
               rprocess=discrete.time.sim(
                 step.fun=ricker.sim
                 ),
               paramnames=c("r","sigma","phi"),
               statenames=c("N","e"),
               measurement.model=y~pois(lambda=N*phi),
               params=c(
                 r=exp(3.8),
                 sigma=0.3,
                 phi=10,
                 N.0=7,
                 e.0=0
                 )
               )
ricker <- simulate(ricker,seed=73691676L)


## ----get-ricker,echo=T,eval=T,results='hide'-----------------------------
pompExample(ricker)


## ----probe-list----------------------------------------------------------
plist <- list(
              probe.marginal("y",ref=obs(ricker),transform=sqrt),
              probe.acf("y",lags=c(0,1,2,3,4),transform=sqrt),
              probe.nlar("y",lags=c(1,1,1,2),powers=c(1,2,3,1),
                         transform=sqrt)
              )


## ----first-probe,eval=F,echo=T-------------------------------------------
## pb.truth <- probe(ricker,probes=plist,nsim=1000,seed=1066L)
## guess <- c(r=20,sigma=1,phi=20,N.0=7,e.0=0)
## pb.guess <- probe(ricker,params=guess,probes=plist,nsim=1000,seed=1066L)

## ----first-probe-eval,eval=T,echo=F--------------------------------------
binary.file <- "ricker-first-probe.rda"
if (file.exists(binary.file)) {
  load(binary.file)
} else {
pb.truth <- probe(ricker,probes=plist,nsim=1000,seed=1066L)
guess <- c(r=20,sigma=1,phi=20,N.0=7,e.0=0)
pb.guess <- probe(ricker,params=guess,probes=plist,nsim=1000,seed=1066L)
  save(pb.truth,pb.guess,guess,file=binary.file,compress='xz')
}


## ----first-probe-plot,eval=F---------------------------------------------
## summary(pb.truth)
## summary(pb.guess)
## plot(pb.truth)
## plot(pb.guess)


## ----ricker-probe-plot,echo=F,results='hide'-----------------------------
binary.file <- "ricker-probe.rda"
if (file.exists(binary.file)) {
  load(binary.file)
} else {
  pb <- probe(ricker,
              probes=list(
                probe.marginal("y",ref=obs(ricker),transform=sqrt),
                probe.acf("y",lags=c(0,1,3),transform=sqrt),
                mean=probe.mean("y",transform=sqrt)
                ),
              transform=TRUE,
              nsim=1000,
              seed=1066L
              )
  save(pb,file=binary.file,compress='xz')
}
plot(pb)


## ----ricker-probe-match-calc,eval=F--------------------------------------
## pm <- probe.match(
##                   pb.guess,
##                   est=c("r","sigma","phi"),
##                   transform=TRUE,
##                   method="Nelder-Mead",
##                   maxit=2000,
##                   seed=1066L,
##                   reltol=1e-8,
##                   trace=3
##                   )
## summary(pm)


## ----ricker-probe.match-eval,echo=F,eval=T,results='hide'----------------
binary.file <- "ricker-probe-match.rda"
if (file.exists(binary.file)) {
  load(binary.file)
} else {
pm <- probe.match(
                  pb.guess,
                  est=c("r","sigma","phi"),
                  transform=TRUE,
                  method="Nelder-Mead",
                  maxit=2000,
                  seed=1066L,
                  reltol=1e-8,
                  trace=3
                  )
summary(pm)
  save(pm,file=binary.file,compress="xz")
}


## ----ricker-mif-calc,eval=F----------------------------------------------
## mf <- mif(
##           ricker,
##           start=guess,
##           Nmif=100,
##           Np=1000,
##           transform=TRUE,
##           cooling.type="geometric",
##           cooling.fraction=0.6,
##           var.factor=2,
##           ic.lag=3,
##           max.fail=50,
##           rw.sd=c(r=0.1,sigma=0.1,phi=0.1)
##           )
## mf <- continue(mf,Nmif=500,max.fail=20)


## ----ricker-mif-eval,echo=F,eval=T,results='hide'------------------------
binary.file <- "ricker-mif.rda"
if (file.exists(binary.file)) {
  load(binary.file)
} else {
mf <- mif(
          ricker,
          start=guess,
          Nmif=100,
          Np=1000,
          transform=TRUE,
          cooling.type="geometric",
          cooling.fraction=0.6,
          var.factor=2,
          ic.lag=3,
          max.fail=50,
          rw.sd=c(r=0.1,sigma=0.1,phi=0.1)
          )
mf <- continue(mf,Nmif=500,max.fail=20)
  save(mf,file=binary.file,compress="xz")
}


## ----ricker-comparison,eval=F,echo=T-------------------------------------
## pf.truth <- pfilter(ricker,Np=1000,max.fail=50,seed=1066L)
## pf.guess <- pfilter(ricker,params=guess,Np=1000,max.fail=50,seed=1066L)
## pf.mf <- pfilter(mf,Np=1000,seed=1066L)
## pf.pm <- pfilter(pm,Np=1000,max.fail=10,seed=1066L)
## pb.mf <- probe(mf,nsim=1000,probes=plist,seed=1066L)
## res <- rbind(
##              cbind(guess=guess,truth=coef(ricker),MLE=coef(mf),PM=coef(pm)),
##              loglik=c(
##                logLik(pf.guess),logLik(pf.truth),logLik(pf.mf),logLik(pf.pm)
##                ),
##              synth.loglik=c(
##                logLik(pb.guess),logLik(pb.truth),logLik(pb.mf),logLik(pm)
##                )
##              )

## ----ricker-comparison-eval,echo=F,eval=T--------------------------------
binary.file <- "ricker-comparison.rda"
if (file.exists(binary.file)) {
  load(binary.file) 
} else {
pf.truth <- pfilter(ricker,Np=1000,max.fail=50,seed=1066L)
pf.guess <- pfilter(ricker,params=guess,Np=1000,max.fail=50,seed=1066L)
pf.mf <- pfilter(mf,Np=1000,seed=1066L)
pf.pm <- pfilter(pm,Np=1000,max.fail=10,seed=1066L)
pb.mf <- probe(mf,nsim=1000,probes=plist,seed=1066L)
res <- rbind(
             cbind(guess=guess,truth=coef(ricker),MLE=coef(mf),PM=coef(pm)),
             loglik=c(
               logLik(pf.guess),logLik(pf.truth),logLik(pf.mf),logLik(pf.pm)
               ),
             synth.loglik=c(
               logLik(pb.guess),logLik(pb.truth),logLik(pb.mf),logLik(pm)
               )
             )
  save(res,file=binary.file,compress='xz')
}

## ----ricker-comparison-show----------------------------------------------
print(res,digits=3)


## ----first-nlf,eval=T,results='hide'-------------------------------------
pompExample(gompertz)
out <- nlf(
           gompertz,
           start=c(r=1,K=2,sigma=0.5,tau=0.5,X.0=1),
           transform.params=TRUE,
           est=c("K","r"),
           lags=c(1,2)
           )


## ----nlf-gompertz-starts,eval=F------------------------------------------
## # pick 5 random starting parameter values
## starts <- replicate(n=5,
##                     {
##                       p <- coef(gompertz)
##                       p[c("K","r")] <- rlnorm(n=2,meanlog=log(p[c("K","r")]),
##                                               sdlog=0.1)
##                       p
##                     },
##                     simplify=FALSE
##                     )


## ----nlf-gompertz-fits,eval=F--------------------------------------------
## out <- list()
## ## Do the fitting.
## ## method, trace, and nasymp are explained below
## for (j in 1:5) {
##   out[[j]] <- nlf(
##                   gompertz,
##                   start=starts[[j]],
##                   transform.data=log,
##                   transform.params=TRUE,
##                   est=c("K","r"),
##                   lags=c(1,2),
##                   seed=7639873L,
##                   method="Nelder-Mead",
##                   trace=4,
##                   skip.se=TRUE,
##                   nasymp=5000
##                   )
## }
## fits <- t(sapply(out,function(x)c(coef(x,c("r","K")),value=logLik(x))))

## ----nlf-fits-eval,echo=F,eval=T,results='hide'--------------------------
binary.file <- "nlf-fits.rda"
if (file.exists(binary.file)) {
  load(binary.file)
} else {
  # pick 5 random starting parameter values
  starts <- replicate(n=5,
                      {
                        p <- coef(gompertz)
                        p[c("K","r")] <- rlnorm(n=2,meanlog=log(p[c("K","r")]),
                                                sdlog=0.1)
                        p
                      },
                      simplify=FALSE
                      )
  out <- list()
  ## Do the fitting. 
  ## method, trace, and nasymp are explained below   
  for (j in 1:5) {
    out[[j]] <- nlf(
                    gompertz,
                    start=starts[[j]],
                    transform.data=log,
                    transform.params=TRUE,
                    est=c("K","r"),
                    lags=c(1,2),
                    seed=7639873L,
                    method="Nelder-Mead",
                    trace=4,
                    skip.se=TRUE,
                    nasymp=5000
                    )  
  }
  fits <- t(sapply(out,function(x)c(coef(x,c("r","K")),value=logLik(x))))
  save(starts,out,fits,file=binary.file,compress="xz")
}


## ----eval=T--------------------------------------------------------------
fits


## ----nlf-my-pomp,eval=T--------------------------------------------------
long.gomp <- simulate(gompertz,times=1:1000)
theta <- coef(long.gomp)


## ----nlf-lag-test-log.r,eval=F-------------------------------------------
## lags <- list(1,2,c(1,2),c(2,3))
## r.vals <- theta["r"]*exp(seq(-0.69,0.69,length=25))
## fvals <- matrix(nrow=25,ncol=4)
## for (j in 1:25) {
##   pars <- theta
##   pars["r"] <- r.vals[j]
##   for(k in 1:4) {
##     fit <- nlf(
##                long.gomp,
##                start=pars,
##                nasymp=5000,
##                lags=lags[[k]],
##                eval.only=TRUE
##                )
##     fvals[j,k] <- logLik(fit)
##   }
## }

## ----nlf-lag-test-log.K,eval=F,echo=F------------------------------------
## K.vals <- theta["K"]*exp(seq(-0.15,0.15,length=25))
## fvals2 <- matrix(nrow=25,ncol=4)
## for (j in 1:25) {
##   pars <- theta
##   pars["K"] <- pars["X.0"] <- K.vals[j]
##   for(k in 1:4) {
##     fit <- nlf(
##                long.gomp,
##                start=pars,
##                nasymp=5000,
##                lags=lags[[k]],
##                eval.only=TRUE
##                )
##     fvals2[j,k] <- logLik(fit)
##   }
## }

## ----nlf-lag-tests,eval=T,echo=F,results='hide'--------------------------
binary.file <- "nlf-lag-tests.rda"
if (file.exists(binary.file)) {
  load(binary.file)
} else {
  long.gomp <- simulate(gompertz,times=1:1000)
  theta <- coef(long.gomp)
  lags <- list(1,2,c(1,2),c(2,3))
  r.vals <- theta["r"]*exp(seq(-0.69,0.69,length=25))
  fvals <- matrix(nrow=25,ncol=4)
  for (j in 1:25) {
    pars <- theta
    pars["r"] <- r.vals[j]
    for(k in 1:4) {
      fit <- nlf(
                 long.gomp,
                 start=pars,
                 nasymp=5000,
                 lags=lags[[k]],
                 eval.only=TRUE
                 )
      fvals[j,k] <- logLik(fit)
    }
  }
  K.vals <- theta["K"]*exp(seq(-0.15,0.15,length=25))
  fvals2 <- matrix(nrow=25,ncol=4)
  for (j in 1:25) {
    pars <- theta
    pars["K"] <- pars["X.0"] <- K.vals[j]
    for(k in 1:4) {
      fit <- nlf(
                 long.gomp,
                 start=pars,
                 nasymp=5000,
                 lags=lags[[k]],
                 eval.only=TRUE
                 )
      fvals2[j,k] <- logLik(fit)
    }
  }
  save(theta,lags,r.vals,K.vals,fvals,fvals2,file=binary.file,compress="xz")
}


## ----nlf-gompertz-plot,fig.height=4,fig.width=6,echo=F-------------------
fvals <- scale(fvals,center=apply(fvals,2,max),scale=FALSE) 
fvals2 <- scale(fvals2,center=apply(fvals2,2,max),scale=FALSE)
op <- par(mfrow=c(1,2),mar=c(5,5,1,1))
matplot(
        r.vals,
        fvals,
        xlab="r",
        lty=1,
        col=c("black","green3","red","purple"),
        ylab="NLF objective function",
        type="o",
        pch=16
        )
abline(v=theta["r"],col="blue")
legend("bottomright",legend=c("1","2","1,2","2,3"),col=c("black","green3","red","purple"),lty=1,pch=16,cex=1,bty="n")

matplot(
        K.vals,
        fvals2,
        xlab="K",
        lty=1,
        col=c("black","green3","red","purple"),
        ylab="NLF objective function",
        type="o",
        pch=16
        )
abline(v=theta["K"],col="blue")
par(op)


## ----nlf-multi-short,eval=F----------------------------------------------
## nreps <- 100
## ndata <- 60
## fvals <- matrix(nrow=nreps,ncol=length(lags))
## new.pomp <- simulate(gompertz,times=1:ndata,nsim=nreps,seed=NULL) # nreps simulated data sets
## for (j in 1:nreps) {
##   for (k in seq_along(lags)) {
##     fit <- nlf(
##                new.pomp[[j]],
##                start=coef(gompertz),
##                nasymp=5000,
##                lags=lags[[k]],
##                eval.only=TRUE
##                )
##     fvals[j,k] <- logLik(fit)
##   }
## }
## fvals <- exp(fvals/ndata)

## ----nlf-multi-short-eval,echo=F,eval=T,results='hide'-------------------
binary.file <- "nlf-multi-short.rda"
if (file.exists(binary.file)) {
  load(binary.file)
} else {
  nreps <- 100
  ndata <- 60
  fvals <- matrix(nrow=nreps,ncol=length(lags))
  new.pomp <- simulate(gompertz,times=1:ndata,nsim=nreps,seed=NULL) # nreps simulated data sets 
  for (j in 1:nreps) {
    for (k in seq_along(lags)) {
      fit <- nlf(
                 new.pomp[[j]], 
                 start=coef(gompertz), 
                 nasymp=5000, 
                 lags=lags[[k]],
                 eval.only=TRUE
                 ) 
      fvals[j,k] <- logLik(fit)
    }
  }
  fvals <- exp(fvals/ndata)
  save(lags,nreps,ndata,fvals,file=binary.file,compress="xz")
}


## ----eval=T--------------------------------------------------------------
apply(fvals,2,function(x)sd(x)/mean(x))


## ----nlf-fit-from-truth,eval=F-------------------------------------------
## true.fit <- nlf(
##                 gompertz,
##                 transform.params=TRUE,
##                 est=c("K","r"),
##                 lags=2,
##                 seed=7639873,
##                 method="Nelder-Mead",
##                 trace=4,
##                 nasymp=5000
##                 )

## ----nlf-fit-from-truth-eval,echo=F,eval=T,results='hide'----------------
binary.file <- "nlf-fit-from-truth.rda"
if (file.exists(binary.file)) {
  load(binary.file)
} else {
  true.fit <- nlf(
                  gompertz,
                  transform.params=TRUE,
                  est=c("K","r"),
                  lags=2,
                  seed=7639873,
                  method="Nelder-Mead",
                  trace=4,
                  nasymp=5000
                  )
  save(true.fit,file=binary.file,compress="xz")
}


## ----echo=F--------------------------------------------------------------
set.seed(32329L)

## ----nlf-boot,eval=F-----------------------------------------------------
## lags <- 2
## ndata <- length(obs(gompertz))
## nboot <- ndata-max(lags)
## nreps <- 100
## pars <- matrix(0,nreps,2)
## bootsamp <- replicate(n=nreps,sample(nboot,replace=TRUE))
## for (j in seq_len(nreps)) {
##   fit <- nlf(
##              gompertz,
##              start=coef(gompertz),
##              transform.params=TRUE,
##              est=c("K","r"),
##              lags=lags,
##              seed=7639873,
##              bootstrap=TRUE,
##              bootsamp=bootsamp[,j],
##              skip.se=TRUE,
##              method="Nelder-Mead",
##              trace=4,
##              nasymp=5000
##              )
##    pars[j,] <- coef(fit,c("r","K"))
## }
## colnames(pars) <- c("r","K")

## ----nlf-boot-eval,echo=F,eval=T,results='hide'--------------------------
binary.file <- "nlf-boot.rda"
if (file.exists(binary.file)) {
  load(binary.file)
} else {
  lags <- 2
  ndata <- length(obs(gompertz))
  nboot <- ndata-max(lags)
  nreps <- 100
  pars <- matrix(0,nreps,2)
  bootsamp <- replicate(n=nreps,sample(nboot,replace=TRUE))
  for (j in seq_len(nreps)) {
    fit <- nlf(
               gompertz,
               start=coef(gompertz),
               transform.params=TRUE,
               est=c("K","r"),
               lags=lags,
               seed=7639873, 
               bootstrap=TRUE, 
               bootsamp=bootsamp[,j],
               skip.se=TRUE, 
               method="Nelder-Mead",
               trace=4,
               nasymp=5000
               )
     pars[j,] <- coef(fit,c("r","K"))
  }
  colnames(pars) <- c("r","K")
  save(pars,file=binary.file,compress="xz")
}

## ------------------------------------------------------------------------
apply(pars,2,sd)


## ----block-bootsamp,eval=F-----------------------------------------------
## bootsamp <- replicate(n=nreps,sample(nboot,size=floor(nboot/3),replace=TRUE))
## bootsamp <- rbind(bootsamp,bootsamp+1,bootsamp+2)


## ----nlf-block-boot,eval=F,echo=F----------------------------------------
## lags <- 2
## ndata <- length(obs(gompertz))
## nboot <- ndata-max(lags)
## nreps <- 100
## pars <- matrix(0,nreps,2)
## bootsamp <- replicate(
##                       n=nreps,
##                       sample(nboot-2,size=floor(nboot/3),replace=TRUE)
##                       )
## bootsamp <- rbind(bootsamp,bootsamp+1,bootsamp+2)
## for (j in seq_len(nreps)) {
##   fit <- nlf(
##              gompertz,
##              transform.params=TRUE,
##              est=c("K","r"),
##              lags=lags,
##              seed=7639873L,
##              bootstrap=TRUE,
##              bootsamp=bootsamp[,j],
##              skip.se=TRUE,
##              method="Nelder-Mead",
##              trace=4,
##              nasymp=5000
##              )
##    pars[j,] <- coef(fit,c("r","K"))
## }
## colnames(pars) <- c("r","K")


## ----nlf-block-boot-eval,eval=T,echo=F,results='hide'--------------------
binary.file <- "nlf-block-boot.rda"
if (file.exists(binary.file)) {
  load(binary.file)
} else {
  lags <- 2
  ndata <- length(obs(gompertz))
  nboot <- ndata-max(lags)
  nreps <- 100
  pars <- matrix(0,nreps,2)
  bootsamp <- replicate(
                        n=nreps,
                        sample(nboot-2,size=floor(nboot/3),replace=TRUE)
                        )
  bootsamp <- rbind(bootsamp,bootsamp+1,bootsamp+2)
  for (j in seq_len(nreps)) {
    fit <- nlf(
               gompertz,
               transform.params=TRUE,
               est=c("K","r"),
               lags=lags,
               seed=7639873L,
               bootstrap=TRUE, 
               bootsamp=bootsamp[,j],
               skip.se=TRUE, 
               method="Nelder-Mead",
               trace=4,
               nasymp=5000
               )
     pars[j,] <- coef(fit,c("r","K"))
  }
  colnames(pars) <- c("r","K")
  save(pars,file=binary.file,compress="xz")
}


## ----bsmc-example-flat-prior-1,echo=T,eval=T,results='hide'--------------
pompExample(ricker)

ricker <- pomp(
               ricker,
               rprior=function (params, ...) {
                 params["r"] <- exp(runif(n=1,min=2,max=5))
                 params["sigma"] <- runif(n=1,min=0.1,max=1)
                 params
               }
               )


## ----bsmc-example-flat-prior-3,eval=F------------------------------------
##   fit1 <- bsmc(ricker,Np=10000,transform=TRUE,
##                est=c("r","sigma"),smooth=0.2,
##                seed=1050180387L)

## ----bsmc-example-flat-prior-eval,eval=T,echo=F--------------------------
binary.file <- "bsmc-ricker-flat-prior.rda"
if (file.exists(binary.file)) {
  load(binary.file)
} else {
pompExample(ricker)

ricker <- pomp(
               ricker,
               rprior=function (params, ...) {
                 params["r"] <- exp(runif(n=1,min=2,max=5))
                 params["sigma"] <- runif(n=1,min=0.1,max=1)
                 params
               }
               )
  fit1 <- bsmc(ricker,Np=10000,transform=TRUE,
               est=c("r","sigma"),smooth=0.2,
               seed=1050180387L)
  save(fit1,file=binary.file,compress="xz")
}


## ----bsmc-example-flat-prior-coef----------------------------------------
signif(coef(fit1),digits=2)


## ----bsmc-example-flat-prior-plot,fig.height=6,fig.width=6,echo=F--------
plot(fit1,pars=c("r","sigma"),thin=5000)


## ----bsmc-example-normal-prior,eval=F,echo=T-----------------------------
## ricker <- pomp(ricker,
##                rprior=function (params, ...) {
##                  x <- rlnorm(n=2,meanlog=c(4,log(0.5)),sdlog=c(3,5))
##                  params[c("r","sigma")] <- x
##                  params
##                }
##                )
## 
## fit2 <- bsmc(ricker,transform=TRUE,Np=10000,
##              est=c("r","sigma"),smooth=0.2,
##              seed=90348704L)
## 

## ----bsmc-example-normal-prior-eval,eval=T,echo=F------------------------
binary.file <- "bsmc-ricker-normal-prior.rda"
if (file.exists(binary.file)) {
  load(binary.file)
} else {
ricker <- pomp(ricker,
               rprior=function (params, ...) {
                 x <- rlnorm(n=2,meanlog=c(4,log(0.5)),sdlog=c(3,5))
                 params[c("r","sigma")] <- x
                 params
               }
               )

fit2 <- bsmc(ricker,transform=TRUE,Np=10000,
             est=c("r","sigma"),smooth=0.2,
             seed=90348704L)

  save(fit2,file=binary.file,compress="xz")
}

## ----bsmc-example-normal-prior-show,echo=T,eval=T------------------------
signif(coef(fit2),digits=2)


## ----sir-step-R----------------------------------------------------------
require(pomp)

sir.step <- function (x, t, params, delta.t, ...) {
  ## unpack the parameters
  N <- params["N"]             # population size
  gamma <- params["gamma"]     # recovery rate
  mu <- params["mu"]           # birth rate = death rate
  beta <- params["beta"]       # contact rate
  foi <- beta*x["I"]/N         # the force of infection
  trans <- c(
             ## births are assumed to be Poisson:
             rpois(n=1,lambda=mu*N*delta.t),
             ## exits from S:
             reulermultinom(n=1,size=x["S"],rate=c(foi,mu),dt=delta.t),
             ## exits from I:
             reulermultinom(n=1,size=x["I"],rate=c(gamma,mu),dt=delta.t),
             ## exits from R:
             reulermultinom(n=1,size=x["R"],rate=c(mu),dt=delta.t)
             )
  ## now connect the compartments
  x["S"] <- x["S"]+trans[1]-trans[2]-trans[3]
  x["I"] <- x["I"]+trans[2]-trans[4]-trans[5]
  x["R"] <- x["R"]+trans[4]-trans[6]
  x["cases"] <- x["cases"]+trans[4]
  x
}


## ----sir-step-C----------------------------------------------------------
sir.step <- '
  double rate[6];		// transition rates
  double trans[6];		// transition numbers

  // compute the transition rates
  rate[0] = mu*N;		// birth into susceptible class
  rate[1] = beta*I/N;           // force of infection
  rate[2] = mu;			// death from susceptible class
  rate[3] = gamma;		// recovery
  rate[4] = mu;			// death from infectious class
  rate[5] = mu; 		// death from recovered class

  // compute the transition numbers
  trans[0] = rpois(rate[0]*dt);	// births are Poisson
  reulermultinom(2,S,&rate[1],dt,&trans[1]);
  reulermultinom(2,I,&rate[3],dt,&trans[3]);
  reulermultinom(1,R,&rate[5],dt,&trans[5]);

  // balance the equations
  S += trans[0]-trans[1]-trans[2];
  I += trans[1]-trans[3]-trans[4];
  R += trans[3]-trans[5];
  incid += trans[3];		// incidence is cumulative recoveries;
'



## ----sir-pomp-def,eval=T,echo=T,results='hide'---------------------------
rmeas <- '
  cases = rnbinom_mu(theta,rho*incid);
'

dmeas <- '
  lik = dnbinom_mu(cases,theta,rho*incid,give_log);
'

pomp(
     data=data.frame(
       cases=NA,
       time=seq(0,10,by=1/52)
       ),
     times="time",
     t0=-1/52,
     dmeasure=Csnippet(dmeas),
     rmeasure=Csnippet(rmeas),
     rprocess=euler.sim(
       step.fun=Csnippet(sir.step),
       delta.t=1/52/20
       ),
     statenames=c("S","I","R","incid"),
     paramnames=c(
       "gamma","mu","theta","beta",
       "N","rho",
       "S.0","I.0","R.0"
       ), 
     zeronames=c("incid"),
     initializer=function(params, t0, ...) {
       x0 <- c(S=0,I=0,R=0,incid=0)
       fracs <- params[c("S.0","I.0","R.0")]
       x0[1:3] <- round(params['N']*fracs/sum(fracs))
       x0
     },
     params=c(
       N=500000,beta=400,
       gamma=26,mu=1/50,rho=0.1,theta=100,
       S.0=26/400,I.0=0.002,R.0=1
       )
     ) -> sir

simulate(sir,seed=1914679908L) -> sir


## ----sir-plot,echo=F-----------------------------------------------------
plot(sir,var=c("cases","incid","S","I"))


## ----seas-basis----------------------------------------------------------
tbasis <- seq(-1,21,by=1/52)
basis <- periodic.bspline.basis(tbasis,nbasis=3,degree=2,period=1,
                                names="seas%d")


## ----complex-sir-def,echo=T,eval=T,results='hide'------------------------
seas.sir.step <- '
  double rate[6];		// transition rates
  double trans[6];		// transition numbers
  double beta;			// transmission rate
  double dW;			// white noise increment
  int k;

  // seasonality in transmission
  for (k = 0, beta = 0.0; k < nbasis; k++)
     beta += (&beta1)[k]*(&seas1)[k];

  // compute the environmental stochasticity
  dW = rgammawn(sigma,dt);

  // compute the transition rates
  rate[0] = mu*N;		// birth into susceptible class
  rate[1] = (iota+beta*I)/N*dW/dt; // force of infection
  rate[2] = mu;			// death from susceptible class
  rate[3] = gamma;		// recovery
  rate[4] = mu;			// death from infectious class
  rate[5] = mu; 		// death from recovered class

  // compute the transition numbers
  trans[0] = rpois(rate[0]*dt);	// births are Poisson
  reulermultinom(2,S,&rate[1],dt,&trans[1]);
  reulermultinom(2,I,&rate[3],dt,&trans[3]);
  reulermultinom(1,R,&rate[5],dt,&trans[5]);

  // balance the equations
  S += trans[0]-trans[1]-trans[2];
  I += trans[1]-trans[3]-trans[4];
  R += trans[3]-trans[5];
  incid += trans[3];	// incidence is cumulative recoveries
'


## ----other-codes,results='hide'------------------------------------------
seas.sir.skel <- '
  double rate[6];		// transition rates
  double term[6];		// transition numbers
  double beta;			// transmission rate
  double dW;			// white noise increment
  int k;
  
  for (k = 0, beta = 0.0; k < nbasis; k++)
     beta += (&beta1)[k]*(&seas1)[k];

  // compute the transition rates
  rate[0] = mu*N;		// birth into susceptible class
  rate[1] = (iota+beta*I)/N;    // force of infection
  rate[2] = mu;			// death from susceptible class
  rate[3] = gamma;		// recovery
  rate[4] = mu;			// death from infectious class
  rate[5] = mu; 		// death from recovered class

  // compute the several terms
  term[0] = rate[0];
  term[1] = rate[1]*S;
  term[2] = rate[2]*S;
  term[3] = rate[3]*I;
  term[4] = rate[4]*I;
  term[5] = rate[5]*R;

  // assemble the differential equations
  DS = term[0]-term[1]-term[2];
  DI = term[1]-term[3]-term[4];
  DR = term[3]-term[5];
  Dincid = term[3]; // accumulate the new I->R transitions
'

partrans <- '
  int k;
  Tgamma = exp(gamma);
  Tmu = exp(mu);
  Tiota = exp(iota);
  for (k = 0; k < nbasis; k++)
    (&Tbeta1)[k] = exp((&beta1)[k]);
  Tsigma = exp(sigma);
  Trho = expit(rho);
  Ttheta = exp(theta);
  from_log_barycentric(&TS_0,&S_0,3);
'

paruntrans <- '
  int k;
  Tgamma = log(gamma);
  Tmu = log(mu);
  Tiota = log(iota);
  for (k = 0; k < nbasis; k++)
    (&Tbeta1)[k] = log((&beta1)[k]);
  Tsigma = log(sigma);
  Trho = logit(rho);
  Ttheta = log(theta);
  to_log_barycentric(&TS_0,&S_0,3);
'

simulate(
         pomp(
              sir,
              rmeasure=Csnippet(rmeas),
              dmeasure=Csnippet(dmeas),
              rprocess=euler.sim(
                step.fun=Csnippet(seas.sir.step),
                delta.t=1/52/20
                ),
              covar=basis,
              tcovar=tbasis,
              skeleton=Csnippet(seas.sir.skel),
              skeleton.type='vectorfield',
              parameter.transform=Csnippet(partrans),
              parameter.inv.transform=Csnippet(paruntrans),
              statenames=c("S","I","R","incid"),
              paramnames=c(
                "gamma","mu","iota","beta1","sigma",
                "N","rho","theta","S.0","I.0","R.0"
                ), 
              zeronames="incid",
              globals="int nbasis = 3;\n",
              initializer=function(params, t0, ...) {
                x0 <- c(S=0,I=0,R=0,incid=0)
                fracs <- params[c("S.0","I.0","R.0")]
                x0[1:3] <- round(params['N']*fracs/sum(fracs))
                x0
              }
              ),
         params=c(
           N=500000,beta1=60,beta2=10,beta3=110,
           gamma=8,mu=1/50,rho=0.5,theta=30,
           iota=20,sigma=0.1,
           S.0=0.13,I.0=0.003,R.0=0.87
           ),
         seed=334849254L
         ) -> complex.sir



## ----seas-basis-plot,echo=F,fig.height=4,fig.width=6---------------------
op <- par(mar=c(5,5,1,5))
matplot(tbasis,basis,xlim=c(0,2),type='l',lwd=2,bty='u',
        lty=1,col=c("red","blue","orange"),xlab="time (yr)",
        ylab=quote("basis functions"~list(s[1],s[2],s[3])))
bb <- coef(complex.sir,c("beta1","beta2","beta3"))
plot.window(c(0,2),range(bb))
lines(tbasis,basis%*%bb,col="black",lwd=3,lty=1)
lines(tbasis,basis%*%bb,col="white",lwd=2.5,lty="23")
axis(side=4)
mtext(
      side=4,
      line=3,
      text=bquote(
        beta==.(beta1)*s[1]+.(beta2)*s[2]+.(beta3)*s[3],
        where=as.list(coef(complex.sir))
        )
      )
par(op)



## ----complex-sir-plot,echo=F---------------------------------------------
plot(complex.sir)


