### R code from vignette source 'intro_to_pomp.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: set-opts
###################################################
  glop <- options(keep.source=TRUE,width=60,continue=" ",prompt=" ")
  library(pomp)
  pdf.options(useDingbats=FALSE)
  set.seed(5384959)


###################################################
### code chunk number 2: gompertz-proc-sim-def
###################################################
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


###################################################
### code chunk number 3: gompertz-meas-sim-def
###################################################
gompertz.meas.sim <- function (x, t, params, ...) {
  ## unpack the parameters:
  tau <- params["tau"]
  ## state at time t:
  X <- x["X"]
  ## generate a simulated observation:
  y <- c(Y=unname(rlnorm(n=1,meanlog=log(X),sd=tau)))
  return(y)
}


###################################################
### code chunk number 4: gompertz-meas-dens-def
###################################################
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


###################################################
### code chunk number 5: first-pomp-construction (eval = FALSE)
###################################################
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


###################################################
### code chunk number 6: set-params
###################################################
theta <- c(r=0.1,K=1,sigma=0.1,tau=0.1,X.0=1)


###################################################
### code chunk number 7: gompertz-first-simulation (eval = FALSE)
###################################################
## gompertz <- simulate(gompertz,params=theta)


###################################################
### code chunk number 8: gompertz-get-data
###################################################
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


###################################################
### code chunk number 9: intro_to_pomp.Rnw:326-327 (eval = FALSE)
###################################################
## plot(gompertz,variables="Y")


###################################################
### code chunk number 10: gompertz-plot
###################################################
plot(gompertz,variables=c("Y"))


###################################################
### code chunk number 11: second-pomp-construction
###################################################
gompertz <- pomp(
                 gompertz,
                 dmeasure=gompertz.meas.dens
                 )


###################################################
### code chunk number 12: gompertz-pfilter-truth (eval = FALSE)
###################################################
## pf <- pfilter(gompertz,params=theta,Np=1000)
## loglik.truth <- logLik(pf)
## loglik.truth


###################################################
### code chunk number 13: gompertz-pfilter-truth-eval
###################################################
set.seed(457645443L)
pf <- pfilter(gompertz,params=theta,Np=1000)
loglik.truth <- logLik(pf)
loglik.truth


###################################################
### code chunk number 14: gompertz-pfilter-truth-alt1 (eval = FALSE)
###################################################
## pf <- pfilter(gompertz,params=coef(gompertz),Np=1000)


###################################################
### code chunk number 15: gompertz-pfilter-truth-alt2 (eval = FALSE)
###################################################
## pf <- pfilter(gompertz,Np=1000)


###################################################
### code chunk number 16: gompertz-pfilter-guess (eval = FALSE)
###################################################
## theta.true <- coef(gompertz)
## theta.guess <- theta.true
## theta.guess[c("r","K","sigma")] <- 1.5*theta.true[c("r","K","sigma")]
## pf <- pfilter(gompertz,params=theta.guess,Np=1000)
## loglik.guess <- logLik(pf)


###################################################
### code chunk number 17: gompertz-pfilter-guess-eval
###################################################
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


###################################################
### code chunk number 18: kalman-filter-def
###################################################
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
    cond.loglik[k] <- dnorm(x=log(Y[k]),mean=M,sd=sqrt(q),log=TRUE)
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


###################################################
### code chunk number 19: kalman-filter-run
###################################################
y <- obs(gompertz)
x0 <- init.state(gompertz)
r <- coef(gompertz,"r")
K <- coef(gompertz,"K")
sigma <- coef(gompertz,"sigma")
tau <- coef(gompertz,"tau")
kf <- kalman.filter(y,x0,r,K,sigma,tau)


###################################################
### code chunk number 20: kalman-likelihood-correction
###################################################
loglik.kalman <- kf$loglik-sum(log(obs(gompertz)))


###################################################
### code chunk number 21: intro_to_pomp.Rnw:465-466 (eval = FALSE)
###################################################
## as(gompertz,"data.frame")


###################################################
### code chunk number 22: intro_to_pomp.Rnw:470-473 (eval = FALSE)
###################################################
## obs(gompertz)
## obs(gompertz,"Y")
## time(gompertz)  


###################################################
### code chunk number 23: intro_to_pomp.Rnw:476-477 (eval = FALSE)
###################################################
## time(gompertz) <- 1:10


###################################################
### code chunk number 24: intro_to_pomp.Rnw:480-482 (eval = FALSE)
###################################################
## timezero(gompertz)
## timezero(gompertz) <- -10


###################################################
### code chunk number 25: intro_to_pomp.Rnw:485-487 (eval = FALSE)
###################################################
## time(gompertz,t0=TRUE)  
## time(gompertz,t0=T) <- seq(from=0,to=10,by=1)


###################################################
### code chunk number 26: intro_to_pomp.Rnw:490-491 (eval = FALSE)
###################################################
## window(gompertz,start=3,end=20)


###################################################
### code chunk number 27: intro_to_pomp.Rnw:495-497 (eval = FALSE)
###################################################
## coef(gompertz)
## coef(gompertz,c("sigma","tau")) <- c(1,0)


###################################################
### code chunk number 28: intro_to_pomp.Rnw:501-503 (eval = FALSE)
###################################################
## states(gompertz)
## states(gompertz,"X")


###################################################
### code chunk number 29: loggomp-pomp-construction
###################################################
gompertz <- pomp(
                 gompertz,
                 parameter.transform=function(params,...){
                   exp(params)
                 },
                 parameter.inv.transform=function(params,...){
                   log(params)
                 }
                 )


###################################################
### code chunk number 30: intro_to_pomp.Rnw:542-543
###################################################
coef(gompertz) <- c(r=0.1,K=1,tau=0.1,sigma=0.1,X.0=1)


###################################################
### code chunk number 31: intro_to_pomp.Rnw:546-547
###################################################
coef(gompertz)


###################################################
### code chunk number 32: intro_to_pomp.Rnw:550-552
###################################################
coef(gompertz,transform=TRUE) <- c(r=log(0.1),K=0,tau=log(0.1),
                sigma=log(0.1),X.0=0)


###################################################
### code chunk number 33: intro_to_pomp.Rnw:555-556
###################################################
coef(gompertz,transform=TRUE)


###################################################
### code chunk number 34: intro_to_pomp.Rnw:559-560
###################################################
coef(gompertz)


###################################################
### code chunk number 35: par-trans-inverse-test
###################################################
# use parameter.inv.transform:
theta <- coef(gompertz,transform=TRUE)  
## theta is on the estimation scale
g2 <- gompertz
## use parameter.transform:
coef(g2,transform=TRUE) <- theta
## compare the internal-scale representations:
identical(coef(gompertz),coef(g2))


###################################################
### code chunk number 36: intro_to_pomp.Rnw:580-581 (eval = FALSE)
###################################################
## demo(gompertz)


###################################################
### code chunk number 37: intro_to_pomp.Rnw:608-611
###################################################
pompExample(gompertz)
theta <- coef(gompertz)
theta.true <- theta


###################################################
### code chunk number 38: gompertz-multi-mif-calc (eval = FALSE)
###################################################
## estpars <- c("r","sigma","tau")
## replicate(
##           n=10,
##           {
##             theta.guess <- theta.true
##             theta.guess[estpars] <- rlnorm(
##                                            n=length(estpars),
##                                            meanlog=log(theta.guess[estpars]),
##                                            sdlog=1
##                                            )
##             mif(
##                 gompertz,
##                 Nmif=100,
##                 start=theta.guess,
##                 transform=TRUE,
##                 pars=estpars,
##                 rw.sd=c(r=0.02,sigma=0.02,tau=0.05),
##                 Np=2000,
##                 var.factor=4,
##                 ic.lag=10,
##                 cooling.type="geometric",
##                 cooling.fraction=0.95,
##                 max.fail=10
##                 )
##           }
##           ) -> mf


###################################################
### code chunk number 39: gompertz-post-mif (eval = FALSE)
###################################################
## theta.true <- coef(gompertz)
## theta.mif <- apply(sapply(mf,coef),1,mean)
## loglik.mif <- replicate(n=10,logLik(pfilter(mf[[1]],params=theta.mif,Np=10000)))
## loglik.mif.est <- logmeanexp(loglik.mif,se=TRUE)
## loglik.true <- replicate(n=10,logLik(pfilter(gompertz,params=theta.true,Np=10000)))
## loglik.true.est <- logmeanexp(loglik.true,se=TRUE)


###################################################
### code chunk number 40: gompertz-multi-mif-eval
###################################################
set.seed(334388458L)
binary.file <- "gompertz-multi-mif.rda"
if (file.exists(binary.file)) {
  load(binary.file)
} else {
  tic <- Sys.time()
estpars <- c("r","sigma","tau")
replicate(
          n=10,
          {
            theta.guess <- theta.true
            theta.guess[estpars] <- rlnorm(
                                           n=length(estpars),
                                           meanlog=log(theta.guess[estpars]),
                                           sdlog=1
                                           )
            mif(
                gompertz,
                Nmif=100,
                start=theta.guess,
                transform=TRUE,
                pars=estpars,
                rw.sd=c(r=0.02,sigma=0.02,tau=0.05),
                Np=2000,
                var.factor=4,
                ic.lag=10,
                cooling.type="geometric",
                cooling.fraction=0.95,
                max.fail=10
                )
          }
          ) -> mf
theta.true <- coef(gompertz)
theta.mif <- apply(sapply(mf,coef),1,mean)
loglik.mif <- replicate(n=10,logLik(pfilter(mf[[1]],params=theta.mif,Np=10000)))
loglik.mif.est <- logmeanexp(loglik.mif,se=TRUE)
loglik.true <- replicate(n=10,logLik(pfilter(gompertz,params=theta.true,Np=10000)))
loglik.true.est <- logmeanexp(loglik.true,se=TRUE)
  toc <- Sys.time()
  etime <- toc-tic
  save(
       mf,estpars,
       theta.mif,theta.true,
       loglik.mif.est,loglik.true.est,
       etime,
       file=binary.file,
       compress="xz"
       )
}
rbind(
      mle=c(signif(theta.mif[estpars],3),loglik=round(loglik.mif.est[1],1),loglik.se=signif(loglik.mif.est[2],2)),
      truth=c(signif(theta.true[estpars],3),loglik=round(loglik.true.est[1],1),loglik.se=signif(loglik.true.est[2],2))
      ) -> results.table


###################################################
### code chunk number 41: intro_to_pomp.Rnw:685-686 (eval = FALSE)
###################################################
## theta.true <- coef(gompertz)
## theta.mif <- apply(sapply(mf,coef),1,mean)
## loglik.mif <- replicate(n=10,logLik(pfilter(mf[[1]],params=theta.mif,Np=10000)))
## loglik.mif.est <- logmeanexp(loglik.mif,se=TRUE)
## loglik.true <- replicate(n=10,logLik(pfilter(gompertz,params=theta.true,Np=10000)))
## loglik.true.est <- logmeanexp(loglik.true,se=TRUE)


###################################################
### code chunk number 42: multi-mif-plot
###################################################
op <- par(mfrow=c(4,1),mar=c(3,3,0,0),mgp=c(2,1,0),bty='l')
loglik <- sapply(mf,function(x)conv.rec(x,"loglik"))
log.r <- sapply(mf,function(x)conv.rec(x,"r"))
log.sigma <- sapply(mf,function(x)conv.rec(x,"sigma"))
log.tau <- sapply(mf,function(x)conv.rec(x,"tau"))
matplot(loglik,type='l',lty=1,xlab="",ylab=expression(log~L),xaxt='n',ylim=max(loglik,na.rm=T)+c(-12,3))
matplot(log.r,type='l',lty=1,xlab="",ylab=expression(log~r),xaxt='n')
matplot(log.sigma,type='l',lty=1,xlab="",ylab=expression(log~sigma),xaxt='n')
matplot(log.tau,type='l',lty=1,xlab="MIF iteration",ylab=expression(log~tau))
par(op)


###################################################
### code chunk number 43: mif-plot
###################################################
op <- par(mfrow=c(4,1),mar=c(3,3,0,0),mgp=c(2,1,0),bty='l')
loglik <- sapply(mf,function(x)conv.rec(x,"loglik"))
log.r <- sapply(mf,function(x)conv.rec(x,"r"))
log.sigma <- sapply(mf,function(x)conv.rec(x,"sigma"))
log.tau <- sapply(mf,function(x)conv.rec(x,"tau"))
matplot(loglik,type='l',lty=1,xlab="",ylab=expression(log~L),xaxt='n',ylim=max(loglik,na.rm=T)+c(-12,3))
matplot(log.r,type='l',lty=1,xlab="",ylab=expression(log~r),xaxt='n')
matplot(log.sigma,type='l',lty=1,xlab="",ylab=expression(log~sigma),xaxt='n')
matplot(log.tau,type='l',lty=1,xlab="MIF iteration",ylab=expression(log~tau))
par(op)


###################################################
### code chunk number 44: first-mif-results-table
###################################################
print(results.table)


###################################################
### code chunk number 45: gompertz-skeleton-def
###################################################
gompertz.skel <- function (x, t, params, ...) {
  r <- params["r"]
  K <- params["K"]
  X <- x["X"]
  S <- exp(-r)
  xnew <- c(X=unname(K^(1-S)*X^S))
  return(xnew)
}


###################################################
### code chunk number 46: new-gompertz-pomp-construction
###################################################
new.gompertz <- pomp(
                     gompertz,
                     skeleton.type="map",
                     skeleton=gompertz.skel
                     )
coef(new.gompertz,"X.0") <- 0.1
coef(new.gompertz,"r") <- 0.1
coef(new.gompertz,"tau") <- 0.05
coef(new.gompertz,"sigma") <- 0.05
new.gompertz <- simulate(new.gompertz,seed=88737400L)


###################################################
### code chunk number 47: gompertz-trajmatch-calc (eval = FALSE)
###################################################
## tm <- traj.match(
##                  new.gompertz,
##                  start=coef(new.gompertz),
##                  transform=TRUE,
##                  est=c("r","K","tau","X.0"),
##                  method="Nelder-Mead",
##                  maxit=1000,
##                  reltol=1e-8
##                  )


###################################################
### code chunk number 48: gompertz-trajmatch-eval
###################################################
binary.file <- "gompertz-trajmatch.rda"
if (file.exists(binary.file)) {
  load(binary.file)
} else {
tm <- traj.match(
                 new.gompertz,
                 start=coef(new.gompertz),
                 transform=TRUE,
                 est=c("r","K","tau","X.0"),
                 method="Nelder-Mead",
                 maxit=1000,
                 reltol=1e-8
                 )
  save(new.gompertz,tm,file=binary.file,compress="xz")
}


###################################################
### code chunk number 49: trajmatch-plot
###################################################
op <- par(mfrow=c(1,1),mar=c(3,3,0,0),mgp=c(2,1,0),bty='l')
plot(time(tm),obs(tm,"Y"),xlab="time",ylab=expression(X,Y),type='o')
lines(time(tm),states(tm,"X"),lwd=2)
par(op)


###################################################
### code chunk number 50: ricker-map-defn
###################################################
ricker.sim <- function (x, t, params, delta.t, ...) {
  e <- rnorm(n=1,mean=0,sd=params["sigma"]) 
  xnew <- c(
            params["r"]*x["N"]*exp(-x["N"]+e),
            e
            )
  names(xnew) <- c("N","e")
  xnew
}


###################################################
### code chunk number 51: ricker-pomp
###################################################
ricker <- pomp(
               data=data.frame(time=seq(0,50,by=1),y=NA),
               times="time",
               t0=0,
               rprocess=discrete.time.sim(
                 step.fun=ricker.sim
                 ),
               measurement.model=y~pois(lambda=N*phi)
               )
coef(ricker) <- c(
                  r=exp(3.8),
                  sigma=0.3,
                  phi=10,
                  N.0=7,
                  e.0=0
                  )
ricker <- simulate(ricker,seed=73691676L)


###################################################
### code chunk number 52: get-ricker
###################################################
pompExample(ricker)


###################################################
### code chunk number 53: probe-list
###################################################
plist <- list(
              probe.marginal("y",ref=obs(ricker),transform=sqrt),
              probe.acf("y",lags=c(0,1,2,3,4),transform=sqrt),
              probe.nlar("y",lags=c(1,1,1,2),powers=c(1,2,3,1),transform=sqrt)
              )


###################################################
### code chunk number 54: first-probe (eval = FALSE)
###################################################
## pb.truth <- probe(ricker,probes=plist,nsim=1000,seed=1066L)
## guess <- c(r=20,sigma=1,phi=20,N.0=7,e.0=0)
## pb.guess <- probe(ricker,params=guess,probes=plist,nsim=1000,seed=1066L)


###################################################
### code chunk number 55: first-probe-eval
###################################################
binary.file <- "ricker-first-probe.rda"
if (file.exists(binary.file)) {
  load(binary.file)
} else {
pb.truth <- probe(ricker,probes=plist,nsim=1000,seed=1066L)
guess <- c(r=20,sigma=1,phi=20,N.0=7,e.0=0)
pb.guess <- probe(ricker,params=guess,probes=plist,nsim=1000,seed=1066L)
  save(pb.truth,pb.guess,guess,file=binary.file,compress='xz')
}


###################################################
### code chunk number 56: first-probe-plot (eval = FALSE)
###################################################
## summary(pb.truth)
## summary(pb.guess)
## plot(pb.truth)
## plot(pb.guess)


###################################################
### code chunk number 57: ricker-probe-plot
###################################################
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


###################################################
### code chunk number 58: ricker-probe-match-calc (eval = FALSE)
###################################################
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


###################################################
### code chunk number 59: ricker-probe.match-eval
###################################################
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


###################################################
### code chunk number 60: ricker-mif-calc (eval = FALSE)
###################################################
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


###################################################
### code chunk number 61: ricker-mif-eval
###################################################
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


###################################################
### code chunk number 62: ricker-comparison (eval = FALSE)
###################################################
## pf.truth <- pfilter(ricker,Np=1000,max.fail=50,seed=1066L)
## pf.guess <- pfilter(ricker,params=guess,Np=1000,max.fail=50,seed=1066L)
## pf.mf <- pfilter(mf,Np=1000,seed=1066L)
## pf.pm <- pfilter(pm,Np=1000,max.fail=10,seed=1066L)
## pb.mf <- probe(mf,nsim=1000,probes=plist,seed=1066L)
## res <- rbind(
##              cbind(guess=guess,truth=coef(ricker),MLE=coef(mf),PM=coef(pm)),
##              loglik=c(
##                pf.guess$loglik,
##                pf.truth$loglik,
##                pf.mf$loglik,
##                pf.pm$loglik
##                ),
##              synth.loglik=c(
##                summary(pb.guess)$synth.loglik,
##                summary(pb.truth)$synth.loglik,
##                summary(pb.mf)$synth.loglik,
##                summary(pm)$synth.loglik
##                )
##              )


###################################################
### code chunk number 63: ricker-comparison-eval
###################################################
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
               pf.guess$loglik,
               pf.truth$loglik,
               pf.mf$loglik,
               pf.pm$loglik
               ),
             synth.loglik=c(
               summary(pb.guess)$synth.loglik,
               summary(pb.truth)$synth.loglik,
               summary(pb.mf)$synth.loglik,
               summary(pm)$synth.loglik
               )
             )
  save(res,file=binary.file,compress='xz')
}


###################################################
### code chunk number 64: ricker-comparison-show
###################################################
print(res,digits=3)


###################################################
### code chunk number 65: first-nlf (eval = FALSE)
###################################################
## pompExample(gompertz)
## out <- nlf(
##            gompertz,
##            start=c(r=1,K=2,sigma=0.5,tau=0.5,X.0=1),
##            partrans=TRUE,
##            est=c("K","r"),
##            lags=c(1,2)
##            )


###################################################
### code chunk number 66: nlf-gompertz-starts (eval = FALSE)
###################################################
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


###################################################
### code chunk number 67: nlf-gompertz-fits (eval = FALSE)
###################################################
## out <- list()
## ## Do the fitting. 
## ## method, trace, and nasymp are explained below   
## for (j in 1:5) {
##   out[[j]] <- nlf(
##                   gompertz,
##                   start=starts[[j]],
##                   transform=log,
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
## fits <- t(sapply(out,function(x)c(x$params[c("r","K")],value=x$value)))


###################################################
### code chunk number 68: nlf-fits-eval
###################################################
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
                  transform=log,
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
fits <- t(sapply(out,function(x)c(x$params[c("r","K")],value=x$value)))
  save(starts,out,fits,file=binary.file,compress="xz")
}


###################################################
### code chunk number 69: intro_to_pomp.Rnw:1120-1121
###################################################
fits


###################################################
### code chunk number 70: nlf-my-pomp
###################################################
long.gomp <- simulate(gompertz,times=1:1000)
theta <- coef(long.gomp)


###################################################
### code chunk number 71: nlf-lag-test-log.r (eval = FALSE)
###################################################
## lags <- list(1,2,c(1,2),c(2,3))
## r.vals <- theta["r"]*exp(seq(-0.69,0.69,length=25))
## fvals <- matrix(nrow=25,ncol=4)
## for (j in 1:25) {
##   pars <- theta
##   pars["r"] <- r.vals[j]
##   for(k in 1:4) {
##     fvals[j,k] <- nlf(
##                       long.gomp,
##                       start=pars,
##                       nasymp=5000,
##                       est=NULL,
##                       lags=lags[[k]],
##                       eval.only=TRUE
##                       )
##   }
## }


###################################################
### code chunk number 72: nlf-lag-test-log.K (eval = FALSE)
###################################################
## K.vals <- theta["K"]*exp(seq(-0.15,0.15,length=25))
## fvals2 <- matrix(nrow=25,ncol=4)
## for (j in 1:25) {
##   pars <- theta
##   pars["K"] <- pars["X.0"] <- K.vals[j]
##   for(k in 1:4) {
##     fvals2[j,k] <- nlf(
##                        long.gomp,
##                        start=pars,
##                        nasymp=5000,
##                        est=NULL,
##                        lags=lags[[k]],
##                        eval.only=TRUE
##                        )
##   }
## }


###################################################
### code chunk number 73: nlf-lag-tests
###################################################
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
    fvals[j,k] <- nlf(
                      long.gomp,
                      start=pars,
                      nasymp=5000,
                      est=NULL,
                      lags=lags[[k]],
                      eval.only=TRUE
                      )
  }
}
K.vals <- theta["K"]*exp(seq(-0.15,0.15,length=25))
fvals2 <- matrix(nrow=25,ncol=4)
for (j in 1:25) {
  pars <- theta
  pars["K"] <- pars["X.0"] <- K.vals[j]
  for(k in 1:4) {
    fvals2[j,k] <- nlf(
                       long.gomp,
                       start=pars,
                       nasymp=5000,
                       est=NULL,
                       lags=lags[[k]],
                       eval.only=TRUE
                       )
  }
}
  save(theta,lags,r.vals,K.vals,fvals,fvals2,file=binary.file,compress="xz")
}


###################################################
### code chunk number 74: nlf-gompertz-plot
###################################################
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


###################################################
### code chunk number 75: nlf-multi-short (eval = FALSE)
###################################################
## nreps <- 100
## ndata <- 60
## fvals <- matrix(nrow=nreps,ncol=length(lags))
## new.pomp <- simulate(gompertz,times=1:ndata,nsim=nreps,seed=NULL) # nreps simulated data sets 
## for (j in 1:nreps) {
##   for (k in seq_along(lags)) {
##     fvals[j,k] <- nlf(
##                       new.pomp[[j]], 
##                       start=coef(gompertz), 
##                       nasymp=5000, 
##                       est=NULL,
##                       lags=lags[[k]],
##                       eval.only=TRUE
##                       ) 
##   }
## }
## fvals <- exp(fvals/ndata)


###################################################
### code chunk number 76: nlf-multi-short-eval
###################################################
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
    fvals[j,k] <- nlf(
                      new.pomp[[j]], 
                      start=coef(gompertz), 
                      nasymp=5000, 
                      est=NULL,
                      lags=lags[[k]],
                      eval.only=TRUE
                      ) 
  }
}
fvals <- exp(fvals/ndata)
  save(lags,nreps,ndata,fvals,file=binary.file,compress="xz")
}


###################################################
### code chunk number 77: intro_to_pomp.Rnw:1266-1267
###################################################
apply(fvals,2,function(x)sd(x)/mean(x))


###################################################
### code chunk number 78: nlf-fit-from-truth (eval = FALSE)
###################################################
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


###################################################
### code chunk number 79: nlf-fit-from-truth-eval
###################################################
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


###################################################
### code chunk number 80: intro_to_pomp.Rnw:1305-1306
###################################################
set.seed(32329L)


###################################################
### code chunk number 81: nlf-boot (eval = FALSE)
###################################################
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
##    pars[j,] <- fit$params[c("r","K")]
## }
## colnames(pars) <- c("r","K")


###################################################
### code chunk number 82: nlf-boot-eval
###################################################
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
   pars[j,] <- fit$params[c("r","K")]
}
colnames(pars) <- c("r","K")
  save(pars,file=binary.file,compress="xz")
}


###################################################
### code chunk number 83: intro_to_pomp.Rnw:1343-1344
###################################################
apply(pars,2,sd)


###################################################
### code chunk number 84: block-bootsamp (eval = FALSE)
###################################################
## bootsamp <- replicate(n=nreps,sample(nboot,size=floor(nboot/3),replace=TRUE))
## bootsamp <- rbind(bootsamp,bootsamp+1,bootsamp+2)


###################################################
### code chunk number 85: nlf-block-boot (eval = FALSE)
###################################################
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
##    pars[j,] <- fit$params[c("r","K")]
## }
## colnames(pars) <- c("r","K")


###################################################
### code chunk number 86: nlf-block-boot-eval
###################################################
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
   pars[j,] <- fit$params[c("r","K")]
}
colnames(pars) <- c("r","K")
  save(pars,file=binary.file,compress="xz")
}


###################################################
### code chunk number 87: bsmc-example-flat-prior-1
###################################################
require(pomp)
pompExample(ricker)


###################################################
### code chunk number 88: bsmc-example-flat-prior-2 (eval = FALSE)
###################################################
## set.seed(136872209L)
## Np <- 10000
## prior1 <- parmat(coef(ricker),nrep=Np)
## prior1["r",] <- exp(runif(n=Np,min=2,max=5))
## prior1["sigma",] <- runif(n=Np,min=0.1,max=1)


###################################################
### code chunk number 89: bsmc-example-flat-prior-3 (eval = FALSE)
###################################################
##   fit1 <- bsmc(ricker,params=prior1,transform=TRUE,
##                est=c("r","sigma"),smooth=0.2)


###################################################
### code chunk number 90: bsmc-example-flat-prior-eval
###################################################
binary.file <- "bsmc-ricker-flat-prior.rda"
if (file.exists(binary.file)) {
  load(binary.file)
} else {
require(pomp)
pompExample(ricker)
set.seed(136872209L)
Np <- 10000
prior1 <- parmat(coef(ricker),nrep=Np)
prior1["r",] <- exp(runif(n=Np,min=2,max=5))
prior1["sigma",] <- runif(n=Np,min=0.1,max=1)
  fit1 <- bsmc(ricker,params=prior1,transform=TRUE,
               est=c("r","sigma"),smooth=0.2)
  save(fit1,file=binary.file,compress="xz")
}


###################################################
### code chunk number 91: bsmc-example-flat-prior-coef
###################################################
signif(coef(fit1),digits=2)


###################################################
### code chunk number 92: bsmc-example-flat-prior-plot
###################################################
  plot(fit1,pars=c("r","sigma"),thin=5000)


###################################################
### code chunk number 93: bsmc-example-normal-prior (eval = FALSE)
###################################################
## set.seed(90348704L)
## Np <- 10000
## prior2 <- parmat(coef(ricker),nrep=Np)
## ## log-normal prior on r
## prior2["r",] <- rlnorm(n=Np,meanlog=4,sdlog=3)
## ## log-normal prior on sigma
## prior2["sigma",] <- rlnorm(n=Np,meanlog=log(0.5),sdlog=5)
## fit2 <- bsmc(ricker,params=prior2,transform=TRUE,
##             est=c("r","sigma"),smooth=0.2)


###################################################
### code chunk number 94: bsmc-example-normal-prior-eval
###################################################
binary.file <- "bsmc-ricker-normal-prior.rda"
if (file.exists(binary.file)) {
  load(binary.file)
} else {
set.seed(90348704L)
Np <- 10000
prior2 <- parmat(coef(ricker),nrep=Np)
## log-normal prior on r
prior2["r",] <- rlnorm(n=Np,meanlog=4,sdlog=3)
## log-normal prior on sigma
prior2["sigma",] <- rlnorm(n=Np,meanlog=log(0.5),sdlog=5)
fit2 <- bsmc(ricker,params=prior2,transform=TRUE,
            est=c("r","sigma"),smooth=0.2)
  save(fit2,file=binary.file,compress="xz")
}


###################################################
### code chunk number 95: bsmc-example-normal-prior-show
###################################################
signif(coef(fit2),digits=2)


###################################################
### code chunk number 96: sir-proc-sim-def
###################################################
sir.proc.sim <- function (x, t, params, delta.t, ...) {
  ## unpack the parameters
  N <- params["N"]             # population size
  gamma <- params["gamma"]     # recovery rate
  mu <- params["mu"]           # birth rate = death rate
  beta <- params["beta"]       # contact rate
  foi <- beta*x["I"]/N                       # the force of infection
  trans <- c(
             rpois(n=1,lambda=mu*N*delta.t), # births are assumed to be Poisson
             reulermultinom(n=1,size=x["S"],rate=c(foi,mu),dt=delta.t), # exits from S
             reulermultinom(n=1,size=x["I"],rate=c(gamma,mu),dt=delta.t), # exits from I
             reulermultinom(n=1,size=x["R"],rate=c(mu),dt=delta.t)        # exits from R
             )
  ## now connect the compartments
  x[c("S","I","R","cases")]+c(
      trans[1]-trans[2]-trans[3],
      trans[2]-trans[4]-trans[5],
      trans[4]-trans[6],
      trans[4]                          # accumulate the recoveries
    )
}


###################################################
### code chunk number 97: sir-pomp-def (eval = FALSE)
###################################################
## simulate(
##          pomp(
##               data=data.frame(
##                 time=seq(1/52,15,by=1/52),
##                 reports=NA
##                 ),
##               times="time",
##               t0=0,
##               rprocess=euler.sim(
##                 step.fun=sir.proc.sim,
##                 delta.t=1/52/20
##                 ),
##               measurement.model=reports~binom(size=cases,prob=rho),
##               initializer=function(params, t0, ic.pars, comp.names, ...){
##                 x0 <- c(S=0,I=0,R=0,cases=0)
##                 N <- params["N"]
##                 fracs <- params[ic.pars]
##                 x0[comp.names] <- round(N*fracs/sum(fracs))
##                 x0
##               },
##               zeronames=c("cases"), # 'cases' is an accumulator variable
##               ic.pars=c("S0","I0","R0"), # initial condition parameters
##               comp.names=c("S","I","R") # names of the compartments
##               ),
##          params=c(
##              N=50000,
##              beta=60,gamma=8,mu=1/50,
##              rho=0.6,
##              S0=8/60,I0=0.002,R0=1-8/60-0.001
##            ),
##          seed=677573454L
##          ) -> sir


###################################################
### code chunk number 98: sir-pomp-def-eval
###################################################
binary.file <- "sir-pomp-def.rda"
if (file.exists(binary.file)) {
  load(binary.file)
} else {
simulate(
         pomp(
              data=data.frame(
                time=seq(1/52,15,by=1/52),
                reports=NA
                ),
              times="time",
              t0=0,
              rprocess=euler.sim(
                step.fun=sir.proc.sim,
                delta.t=1/52/20
                ),
              measurement.model=reports~binom(size=cases,prob=rho),
              initializer=function(params, t0, ic.pars, comp.names, ...){
                x0 <- c(S=0,I=0,R=0,cases=0)
                N <- params["N"]
                fracs <- params[ic.pars]
                x0[comp.names] <- round(N*fracs/sum(fracs))
                x0
              },
              zeronames=c("cases"), # 'cases' is an accumulator variable
              ic.pars=c("S0","I0","R0"), # initial condition parameters
              comp.names=c("S","I","R") # names of the compartments
              ),
         params=c(
             N=50000,
             beta=60,gamma=8,mu=1/50,
             rho=0.6,
             S0=8/60,I0=0.002,R0=1-8/60-0.001
           ),
         seed=677573454L
         ) -> sir
  save(sir,file=binary.file,compress='xz')
}


###################################################
### code chunk number 99: sir-pomp-def-with-skel (eval = FALSE)
###################################################
## sir <- pomp(
##             sir,
##               skeleton.type="vectorfield",
##               skeleton=function(t,x,params,...){
##                 N <- params["N"]             # population size
##                 gamma <- params["gamma"]     # recovery rate
##                 mu <- params["mu"]           # birth rate = death rate
##                 beta <- params["beta"]       # contact rate
##                 foi <- beta*x["I"]/N # the force of infection
##                 terms <- c(
##                            mu*N,               # births
##                            x["S"]*c(foi,mu),   # exits from S
##                            x["I"]*c(gamma,mu), # exits from I
##                            x["R"]*mu           # exits from R
##                            )
##                 ## now put the equations together:
##                 c(
##                   cases=unname(terms[4]),               # d(cases)/dt
##                   S=unname(terms[1]-terms[2]-terms[3]), # dS/dt
##                   I=unname(terms[2]-terms[4]-terms[5]), # dI/dt
##                   R=unname(terms[4]-terms[6])           # dR/dt
##                   )
##               }
## )


###################################################
### code chunk number 100: sir-plot
###################################################
plot(sir)


###################################################
### code chunk number 101: seas-basis
###################################################
tbasis <- seq(0,20,by=1/52)
basis <- periodic.bspline.basis(tbasis,nbasis=3,degree=2,period=1,names="seas%d")


###################################################
### code chunk number 102: complex-sir-def (eval = FALSE)
###################################################
## complex.sir.proc.sim <- function (x, t, params, delta.t, covars, ...) {
##   ## unpack the parameters
##   N <- params["N"]                 # population size
##   gamma <- params["gamma"]         # recovery rate
##   mu <- params["mu"]               # birth rate = death rate
##   iota <- params["iota"]           # import rate
##   b <- params[c("b1","b2","b3")]   # contact-rate coefficients
##   beta <- b%*%covars               # flexible seasonality
##   beta.sd <- params["beta.sd"]     # extrademographic noise intensity
##   dW <- rgammawn(n=1,sigma=beta.sd,dt=delta.t) # Gamma white noise
##   foi <- (beta*x["I"]/N+iota)*dW/delta.t # the force of infection
##   trans <- c(
##              rpois(n=1,lambda=mu*N*delta.t), # births are assumed to be Poisson
##              reulermultinom(n=1,size=x["S"],rate=c(foi,mu),dt=delta.t), # exits from S
##              reulermultinom(n=1,size=x["I"],rate=c(gamma,mu),dt=delta.t), # exits from I
##              reulermultinom(n=1,size=x["R"],rate=c(mu),dt=delta.t)        # exits from R
##              )
##   ## now connect the compartments
##   x[c("S","I","R","cases","W")]+
##     c(
##       trans[1]-trans[2]-trans[3],
##       trans[2]-trans[4]-trans[5],
##       trans[4]-trans[6],
##       trans[4],                   # accumulate the recoveries
##       (dW-delta.t)/beta.sd        # mean = 0, var = delta.t
##     )
## }
## 
## simulate(
##          pomp(
##               sir,
##               tcovar=tbasis,
##               covar=basis,
##               rprocess=euler.sim(
##                 complex.sir.proc.sim,
##                 delta.t=1/52/20
##                 ),
##               initializer=function(params, t0, ic.pars, comp.names, ...){
##                 x0 <- c(S=0,I=0,R=0,cases=0,W=0)
##                 N <- params["N"]
##                 fracs <- params[ic.pars]
##                 x0[comp.names] <- round(N*fracs/sum(fracs))
##                 x0
##               }
##               ),
##          params=c(
##              N=50000,
##              b1=60,b2=10,b3=110,
##              gamma=8,mu=1/50,
##              rho=0.6,
##              iota=0.01,beta.sd=0.1,
##              S0=8/60,I0=0.002,R0=1-8/60-0.001
##            ),
##          seed=8274355L
##          ) -> complex.sir


###################################################
### code chunk number 103: complex-sir-def-eval
###################################################
binary.file <- "complex-sir-def.rda"
if (file.exists(binary.file)) {
  load(binary.file)
} else {
complex.sir.proc.sim <- function (x, t, params, delta.t, covars, ...) {
  ## unpack the parameters
  N <- params["N"]                 # population size
  gamma <- params["gamma"]         # recovery rate
  mu <- params["mu"]               # birth rate = death rate
  iota <- params["iota"]           # import rate
  b <- params[c("b1","b2","b3")]   # contact-rate coefficients
  beta <- b%*%covars               # flexible seasonality
  beta.sd <- params["beta.sd"]     # extrademographic noise intensity
  dW <- rgammawn(n=1,sigma=beta.sd,dt=delta.t) # Gamma white noise
  foi <- (beta*x["I"]/N+iota)*dW/delta.t # the force of infection
  trans <- c(
             rpois(n=1,lambda=mu*N*delta.t), # births are assumed to be Poisson
             reulermultinom(n=1,size=x["S"],rate=c(foi,mu),dt=delta.t), # exits from S
             reulermultinom(n=1,size=x["I"],rate=c(gamma,mu),dt=delta.t), # exits from I
             reulermultinom(n=1,size=x["R"],rate=c(mu),dt=delta.t)        # exits from R
             )
  ## now connect the compartments
  x[c("S","I","R","cases","W")]+
    c(
      trans[1]-trans[2]-trans[3],
      trans[2]-trans[4]-trans[5],
      trans[4]-trans[6],
      trans[4],                   # accumulate the recoveries
      (dW-delta.t)/beta.sd        # mean = 0, var = delta.t
    )
}

simulate(
         pomp(
              sir,
              tcovar=tbasis,
              covar=basis,
              rprocess=euler.sim(
                complex.sir.proc.sim,
                delta.t=1/52/20
                ),
              initializer=function(params, t0, ic.pars, comp.names, ...){
                x0 <- c(S=0,I=0,R=0,cases=0,W=0)
                N <- params["N"]
                fracs <- params[ic.pars]
                x0[comp.names] <- round(N*fracs/sum(fracs))
                x0
              }
              ),
         params=c(
             N=50000,
             b1=60,b2=10,b3=110,
             gamma=8,mu=1/50,
             rho=0.6,
             iota=0.01,beta.sd=0.1,
             S0=8/60,I0=0.002,R0=1-8/60-0.001
           ),
         seed=8274355L
         ) -> complex.sir
  save(complex.sir,file=binary.file,compress='xz')
}


###################################################
### code chunk number 104: seas-basis-plot
###################################################
op <- par(mar=c(5,5,1,5))
matplot(tbasis,basis,xlim=c(0,2),type='l',lwd=2,bty='u',
        lty=1,col=c("red","blue","orange"),xlab="time (yr)",
        ylab=quote("basis functions"~list(s[1],s[2],s[3])))
bb <- coef(complex.sir,c("b1","b2","b3"))
plot.window(c(0,2),c(0,1)*max(bb))
lines(tbasis,basis%*%bb,col="black",lwd=3,lty=1)
lines(tbasis,basis%*%bb,col="white",lwd=2.5,lty="23")
axis(side=4)
mtext(
      side=4,
      line=2,
      text=bquote(
        beta==.(b1)*s[1]+.(b2)*s[2]+.(b3)*s[3],
        where=as.list(coef(complex.sir))
        )
      )
par(op)


###################################################
### code chunk number 105: complex-sir-plot
###################################################
plot(complex.sir)


###################################################
### code chunk number 106: restore-opts
###################################################
options(glop)


