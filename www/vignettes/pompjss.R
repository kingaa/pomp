## ----packages,include=F,cache=F------------------------------------------
require(knitr)
require(pomp)
require(coda)
require(foreach)
require(doMC)
require(ggplot2)
require(grid)
require(magrittr)
require(plyr)
require(reshape2)
require(xtable)

stopifnot(packageVersion("pomp")>="0.65-1")


## ----set-opts,include=F,cache=F------------------------------------------
options(
        scipen=2,
        help_type="html",
        stringsAsFactors=FALSE,
        prompt="R> ",
        continue="+  ",
        width=70,
        useFancyQuotes=FALSE,
        knitr.package.concordance=TRUE,
        reindent.spaces=2,
        xtable.comment=FALSE
        )

opts_chunk$set(
               progress=TRUE,prompt=TRUE,highlight=FALSE,
               tidy=TRUE,
               tidy.opts=list(
                 keep.blank.line=FALSE
                 ),
               comment="",
               warning=FALSE,message=FALSE,error=TRUE,
               echo=TRUE,
               strip.white=TRUE,
               results="markup",
               indent=NA,
               background="#FFFFFF00",
               size="normalsize",
               fig.path="figure/pompjss-",
               fig.lp="fig:",
               fig.align="left",
               fig.show="asis",
               fig.height=5,fig.width=8,
               dpi=300,
               dev="png",
               dev.args=list(
                 bg="transparent",
                 pointsize=12
                 )
               )


## ----set-seed,cache=F,include=F------------------------------------------
set.seed(5384959L)
options(cores=10)


## ----bake,include=FALSE,cache=FALSE--------------------------------------
bake <- function (file, expr) {
  expr <- substitute(expr)
  if (file.exists(file)) {
    load(file,envir=parent.frame())
  } else {
    env <- new.env()
    eval(expr,envir=env)
    save(list=ls(envir=env),envir=env,file=file)
    for (n in ls(envir=env))
      assign(n,get(n,envir=env),envir=parent.frame())
  }
  invisible(NULL)
}

## ----timing1,echo=F,cache=F----------------------------------------------
bigtick <- Sys.time()

## ----eval=F--------------------------------------------------------------
## library("pomp")
## pompExample(gompertz)

## ----eval=F--------------------------------------------------------------
## obs(gompertz)
## states(gompertz)
## as.data.frame(gompertz)
## plot(gompertz)
## timezero(gompertz)
## time(gompertz)
## coef(gompertz)
## init.state(gompertz)

## ----gomp1---------------------------------------------------------------
gompertz.proc.sim <- function (x, t, params, delta.t, ...) {
   eps <- exp(rnorm(n=1,mean=0,sd=params["sigma"]))
   S <- exp(-params["r"]*delta.t)
   setNames(params["K"]^(1-S)*x["X"]^S*eps,"X")
 }

## ----gomp2,tidy.opts=list(width.cutoff=72)-------------------------------
gompertz.meas.sim <- function (x, t, params, ...) {
   setNames(rlnorm(n=1,meanlog=log(x["X"]),sd=params["tau"]),"Y")
 }

## ----gomp3,tidy.opts=list(width.cutoff=65)-------------------------------
gompertz.meas.dens <- function (y, x, t, params, log, ...) {
   dlnorm(x=y["Y"],meanlog=log(x["X"]),sdlog=params["tau"],log=log)
 }

## ----gomp5---------------------------------------------------------------
gompertz <- pomp(data=data.frame(time=1:100, Y=NA), times="time",
                 rprocess=discrete.time.sim(step.fun=gompertz.proc.sim,
                   delta.t=1), rmeasure=gompertz.meas.sim, t0=0)

## ----gomp6---------------------------------------------------------------
theta <- c(r=0.1,K=1,sigma=0.1,tau=0.1,X.0=1)

## ----gomp7-setup,echo=F,results="hide"-----------------------------------
set.seed(340398091L)

## ----gomp7---------------------------------------------------------------
gompertz <- simulate(gompertz,params=theta)

## ----gomp8,eval=F--------------------------------------------------------
## plot(gompertz, variables = "Y")

## ----gompertz-plot,echo=F,fig.height=3,fig.width=5-----------------------
op <- par(mar=c(3,3,2,0),mgp=c(2,1,0))
plot(Y~time,data=as.data.frame(gompertz),type="l")
title("gompertz",line=1,cex.main=1)
par(op)

## ----gomp9---------------------------------------------------------------
gompertz <- pomp(gompertz,dmeasure=gompertz.meas.dens)

## ----pfilter1-setup,eval=T,echo=F,results="hide"-------------------------
set.seed(334388458L)

## ----pfilter1-calc,eval=T,cache=T,results="markup",echo=T----------------
pf <- pfilter(gompertz,params=theta,Np=1000)
loglik.truth <- logLik(pf)
loglik.truth

## ----pfilter1-followup,echo=F,results="hide",eval=T,cache=T--------------
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

kalman <- function (x, object, params) {
  # evaluates gompertz likelihood
  # parameters X0 and K are fixed at 1
  # other parameters are taken from x (or, if not in x, from params)
  Y <- obs(object)
  p <- params
  p[names(x)] <- x
  X0 <- 1
  r <- p["r"]
  K <- 1
  sigma <- p["sigma"]
  tau <- p["tau"]
  -kalman.filter(Y, X0, r, K, sigma, tau)$loglik
}

exact.loglik.truth <- -kalman(coef(gompertz),gompertz,coef(gompertz))

## ----pfilter2-calc,eval=FALSE--------------------------------------------
##  pf <- pfilter(gompertz,params=coef(gompertz),Np=1000)

## ----pfilter3-calc,eval=FALSE--------------------------------------------
##  pf <- pfilter(gompertz,Np=1000)

## ----pfilter4-setup,eval=T,echo=F,results="hide"-------------------------
set.seed(334388458L)

## ----pfilter4-calc,eval=T,results="markup",cache=T-----------------------
theta.guess <- theta.true <- coef(gompertz)
theta.guess[c("r","K","sigma")] <- 1.5*theta.true[c("r","K","sigma")]
pf <- pfilter(gompertz,params=theta.guess,Np=1000)
loglik.guess <- logLik(pf)
loglik.guess

## ----gomp4---------------------------------------------------------------
gompertz.log.tf <- function (params, ...) log(params)
gompertz.exp.tf <- function (params, ...) exp(params)

## ----gompertz-transforms,tidy=F------------------------------------------
gompertz <- pomp(gompertz, toEstimationScale = gompertz.log.tf,
                 fromEstimationScale = gompertz.exp.tf)

## ----gompertz-mif-setup,echo=F,results="hide"----------------------------
dat1 <- as.data.frame(gompertz)
pompExample(gompertz)
gompertz <- simulate(window(gompertz,start=1),seed=340398091L)
dat2 <- as.data.frame(gompertz)
stopifnot(all.equal(dat1[c("time","Y")],dat2[c("time","Y")]))
theta <- coef(gompertz)
theta.true <- theta

## ----gompertz-mif-inner1,echo=F,eval=F-----------------------------------
## theta.guess <- theta.true
## rlnorm(n = length(estpars), meanlog = log(theta.guess[estpars]),
##        sdlog = 1) -> theta.guess[estpars]

## ----gompertz-mif-inner2,echo=F,eval=F-----------------------------------
## mif(gompertz, Nmif = 100, start = theta.guess, transform = TRUE,
##     Np = 2000, var.factor = 2, cooling.fraction = 0.7,
##     rw.sd = c(r = 0.02, sigma = 0.02, tau = 0.02))

## ----gompertz-mif-inner3,echo=F,eval=F-----------------------------------
## pf <- replicate(n = 10, logLik(pfilter(mf, Np = 10000)))
## logmeanexp(pf)

## ----gompertz-mif-estpar,eval=F,echo=F-----------------------------------
## estpars <- c("r", "sigma", "tau")

## ----gompertz-mif-post,eval=F,echo=F-------------------------------------
## mf1 <- mif1[[which.max(pf1)]]
## theta.mif <- coef(mf1)
## loglik.mif <- replicate(n = 10, logLik(pfilter(mf1,Np = 10000)))
## loglik.mif <- logmeanexp(loglik.mif,se=TRUE)
## theta.true <- coef(gompertz)
## loglik.true <- replicate(n = 10, logLik(pfilter(gompertz, Np = 20000)))
## loglik.true <- logmeanexp(loglik.true,se=TRUE)

## ----gompertz-mif-demo,eval=F,echo=T,tidy=F------------------------------
## estpars <- c("r", "sigma", "tau")
## library("foreach")
## mif1 <- foreach(i = 1:10, .combine = c) %dopar% {
##   theta.guess <- theta.true
##   rlnorm(n = length(estpars), meanlog = log(theta.guess[estpars]),
##          sdlog = 1) -> theta.guess[estpars]
##   mif(gompertz, Nmif = 100, start = theta.guess, transform = TRUE,
##       Np = 2000, var.factor = 2, cooling.fraction = 0.7,
##       rw.sd = c(r = 0.02, sigma = 0.02, tau = 0.02))
## }
## pf1 <- foreach(mf = mif1, .combine = c) %dopar% {
##   pf <- replicate(n = 10, logLik(pfilter(mf, Np = 10000)))
##   logmeanexp(pf)
## }

## ----gompertz-mif-eval,echo=F,results="hide",cache=F---------------------
bake("gompertz-mif.rda",{

  require(doMC)
  require(foreach)
  registerDoMC()

  save.seed <- .Random.seed
  set.seed(334388458L,kind="L'Ecuyer")

  estpars <- c("r", "sigma", "tau")

  tic <- Sys.time()
  mif1 <- foreach(i=1:10,
    .inorder=FALSE,.packages="pomp",.combine=c,
    .options.multicore=list(set.seed=TRUE)) %dopar% {
    theta.guess <- theta.true
    rlnorm(n = length(estpars), meanlog = log(theta.guess[estpars]),
           sdlog = 1) -> theta.guess[estpars]
    mif(gompertz, Nmif = 100, start = theta.guess, transform = TRUE,
        Np = 2000, var.factor = 2, cooling.fraction = 0.7,
        rw.sd = c(r = 0.02, sigma = 0.02, tau = 0.02))
  }

  pf1 <- foreach(mf=mif1,
  .inorder=TRUE,.packages="pomp",.combine=c,
  .options.multicore=list(set.seed=TRUE)) %dopar% {
    pf <- replicate(n = 10, logLik(pfilter(mf, Np = 10000)))
    logmeanexp(pf)
   }
  toc <- Sys.time()
  mifTime <- toc-tic
  .Random.seed <<- save.seed

  mf1 <- mif1[[which.max(pf1)]]
  theta.mif <- coef(mf1)
  loglik.mif <- replicate(n = 10, logLik(pfilter(mf1,Np = 10000)))
  loglik.mif <- logmeanexp(loglik.mif,se=TRUE)
  theta.true <- coef(gompertz)
  loglik.true <- replicate(n = 10, logLik(pfilter(gompertz, Np = 20000)))
  loglik.true <- logmeanexp(loglik.true,se=TRUE)

  kalm.fit1 <- optim(
    par=theta.guess[estpars],
    fn=kalman,
    object=gompertz,
    params=coef(gompertz),
    hessian=TRUE,
    control=list(trace=2)
  )

  theta.mle <- kalm.fit1$par
  exact.loglik.maximized <- -kalm.fit1$value
  exact.loglik.mif1 <- -kalman(coef(mf1),gompertz,coef(gompertz))

  mle.po <- gompertz
  coef(mle.po,names(theta.mle)) <- unname(theta.mle)
  loglik.mle <- replicate(n=10,logLik(pfilter(mle.po,Np=20000)))
  loglik.mle <- logmeanexp(loglik.mle,se=TRUE)
})

## ----gompertz-mif-results,echo=F,eval=T,results="hide"-------------------
rbind(
      `truth`=c(signif(theta.true[estpars],3),round(loglik.true,2),round(exact.loglik.truth,2)),
      `\\code{mif} MLE`=c(signif(theta.mif[estpars],3),round(loglik.mif,2),round(exact.loglik.mif1,2)),
      `exact MLE`=c(signif(theta.mle[estpars],3),round(loglik.mle,2),round(exact.loglik.maximized,2))
     ) -> results.table
pretty.pars <- c(r="$r$",sigma="$\\sigma$",tau="$\\tau$")
colnames(results.table) <- c(pretty.pars[estpars],"$\\loglikMC$","s.e.","$\\loglik$")

## ----mif3,echo=T,eval=F--------------------------------------------------
## mf1 <- mif1[[which.max(pf1)]]
## theta.mif <- coef(mf1)
## loglik.mif <- replicate(n = 10, logLik(pfilter(mf1,Np = 10000)))
## loglik.mif <- logmeanexp(loglik.mif,se=TRUE)
## theta.true <- coef(gompertz)
## loglik.true <- replicate(n = 10, logLik(pfilter(gompertz, Np = 20000)))
## loglik.true <- logmeanexp(loglik.true,se=TRUE)

## ----mif-plot,echo=F,cache=TRUE,fig.height=6-----------------------------
op <- par(mfrow=c(4,1),mar=c(3,4,0.3,0),mgp=c(2,1,0),
          bty="l",cex.axis=1.2,cex.lab=1.4)
loglik <- do.call(cbind,conv.rec(mif1,"loglik"))
log.r <- do.call(cbind,conv.rec(mif1,"r"))
log.sigma <- do.call(cbind,conv.rec(mif1,"sigma"))
log.tau <- do.call(cbind,conv.rec(mif1,"tau"))
matplot(loglik,type="l",lty=1,xlab="",ylab=expression(log~L),xaxt="n",ylim=max(loglik,na.rm=T)+c(-30,3))
matplot(log.r,type="l",lty=1,xlab="",ylab=expression(log~r),xaxt="n")
matplot(log.sigma,type="l",lty=1,xlab="",ylab=expression(log~sigma),xaxt="n")
matplot(log.tau,type="l",lty=1,xlab="mif iteration",ylab=expression(log~tau))
par(op)

## ----gompertz-multi-mif-table,echo=F,results="asis"----------------------
require(xtable)
options(
xtable.sanitize.text.function=function(x)x,
xtable.floating=FALSE
)
print(xtable(results.table,align="r|cccccc",digits=c(0,4,4,4,2,2,2)))

## ----gompertz-dprior1,tidy=F---------------------------------------------
hyperparams <- list(min = coef(gompertz)/10, max = coef(gompertz) * 10)

## ----gompertz-dprior2,tidy=FALSE-----------------------------------------
gompertz.dprior <- function (params, ..., log) {
  f <- sum(dunif(params, min = hyperparams$min, max = hyperparams$max,
                 log = TRUE))
  if (log) f else exp(f)
}

## ----pmcmc-inner,eval=F,echo=F-------------------------------------------
## pmcmc(pomp(gompertz, dprior = gompertz.dprior), start = theta.mif,
##       Nmcmc = 40000, Np = 100, max.fail = Inf,
##       proposal=mvn.diag.rw(c(r = 0.01, sigma = 0.01, tau = 0.01)))

## ----pmcmc-demo,eval=FALSE,echo=T,tidy=F---------------------------------
## pmcmc1 <- foreach(i=1:5,.combine=c) %dopar% {
##   pmcmc(pomp(gompertz, dprior = gompertz.dprior), start = theta.mif,
##         Nmcmc = 40000, Np = 100, max.fail = Inf,
##         proposal=mvn.diag.rw(c(r = 0.01, sigma = 0.01, tau = 0.01)))
## }

## ----pmcmc-eval,echo=F,results="hide",cache=F----------------------------
require(pomp)
require(coda)

bake("pmcmc.rda",{

  pompExample(gompertz)

  tic <- Sys.time()
  require(doMC)
  require(foreach)
  registerDoMC()

  save.seed <- .Random.seed
  set.seed(334388458L,kind="L'Ecuyer")

  pmcmc1 <- foreach(
    i=1:5,
    .inorder=FALSE,
    .packages="pomp",
    .combine=c,
    .options.multicore=list(set.seed=TRUE)
  ) %dopar%
  {
    pmcmc(pomp(gompertz, dprior = gompertz.dprior), start = theta.mif,
          Nmcmc = 40000, Np = 100, max.fail = Inf,
          proposal=mvn.diag.rw(c(r = 0.01, sigma = 0.01, tau = 0.01)))
  }

toc <- Sys.time()
pmcmcTime <- toc-tic
.Random.seed <<- save.seed

pmcmc.traces <- conv.rec(pmcmc1,c("r","sigma","tau"))
pmcmc.traces <- window(pmcmc.traces,start=20001,thin=40)
ess.pmcmc <- effectiveSize(pmcmc.traces)
rm(pmcmc1,save.seed,tic,toc)
})

## ----pmcmc-diagnostics,results="hide",fig.show="hide",echo=F,eval=T------
gelman.diag(pmcmc.traces)
gelman.plot(pmcmc.traces)
autocorr.plot(pmcmc.traces[[1]])
hist(rle(unlist(pmcmc.traces[,"r"]))$length)

## ----pmcmc-plot,echo=F,eval=T,results="hide",cache=TRUE------------------
op <- par(mar=c(4,3.5,0,1),mfcol=c(3,2),mgp=c(2.5,1,0),cex.axis=1.5,cex.lab=2)
traceplot(pmcmc.traces[,"r"],smooth=TRUE,xlab="",ylab=expression(r),lty=1)
traceplot(pmcmc.traces[,"sigma"],smooth=TRUE,xlab="",ylab=expression(sigma),lty=1)
traceplot(pmcmc.traces[,"tau"],smooth=TRUE,xlab="PMCMC iteration",ylab=expression(tau),lty=1)
densplot(pmcmc.traces[,"r"],show.obs=FALSE,xlab="",main="")
mtext(side=1,line=2,text=expression(r))
abline(v=coef(gompertz,"r"))
densplot(pmcmc.traces[,"sigma"],show.obs=FALSE,xlab="",main="")
mtext(side=1,line=2,text=expression(sigma))
abline(v=coef(gompertz,"sigma"))
densplot(pmcmc.traces[,"tau"],show.obs=FALSE,xlab=expression(tau),main="")
abline(v=coef(gompertz,"tau"))
par(op)

## ----ricker-map-defn,tidy=F----------------------------------------------
ricker.sim <- "
   e = rnorm(0, sigma);
   N = r * N * exp(-N + e);
"
ricker.rmeas <- "
   y = rpois(phi * N);
"
ricker.dmeas <- "
   lik = dpois(y, phi * N, give_log);
"

## ----ricker-trans,tidy=F-------------------------------------------------
log.trans <- "
   Tr = log(r);
   Tsigma = log(sigma);
   Tphi = log(phi);
   TN_0 = log(N_0);"
exp.trans <- "
   Tr = exp(r);
   Tsigma = exp(sigma);
   Tphi = exp(phi);
   TN_0 = exp(N_0);"

## ----ricker-pomp,tidy=F--------------------------------------------------
pomp(data = data.frame(time = seq(0, 50, by = 1), y = NA),
     rprocess = discrete.time.sim(step.fun = Csnippet(ricker.sim),
       delta.t = 1), rmeasure = Csnippet(ricker.rmeas),
     dmeasure = Csnippet(ricker.dmeas),
     toEstimationScale = Csnippet(log.trans),
     fromEstimationScale = Csnippet(exp.trans),
     paramnames = c("r", "sigma", "phi", "N.0", "e.0"),
     statenames = c("N", "e"), times = "time", t0 = 0,
     params = c(r = exp(3.8), sigma = 0.3, phi = 10,
       N.0 = 7, e.0 = 0)) -> ricker
ricker <- simulate(ricker,seed=73691676L)

## ----get-ricker,echo=F,eval=T,results="hide"-----------------------------
dat1 <- as.data.frame(ricker)
pompExample(ricker)
dat2 <- as.data.frame(ricker)
stopifnot(all.equal(dat1[c("time","y")],dat2[c("time","y")]))

## ----probe-list,tidy=FALSE-----------------------------------------------
plist <- list(probe.marginal("y", ref = obs(ricker), transform = sqrt),
              probe.acf("y", lags = c(0, 1, 2, 3, 4), transform = sqrt),
              probe.nlar("y", lags = c(1, 1, 1, 2), powers = c(1, 2, 3, 1),
                         transform = sqrt))

## ----first-probe,eval=T,echo=T,cache=T-----------------------------------
pb.truth <- probe(ricker,probes=plist,nsim=1000,seed=1066L)
guess <- c(r=20,sigma=1,phi=20,N.0=7,e.0=0)
pb.guess <- probe(ricker,params=guess,probes=plist,nsim=1000,seed=1066L)

## ----first-probe-plot,eval=F---------------------------------------------
## summary(pb.truth)
## summary(pb.guess)
## plot(pb.truth)
## plot(pb.guess)

## ----ricker-probe-plot,echo=F,cache=T,results="hide",dpi=600,dev.args=list(bg="transparent",pointsize=9),fig.height=4,fig.width=4,out.width="\\textwidth"----
  pb <- probe(ricker,
              probes=list(
                probe.marginal("y",ref=obs(ricker),transform=sqrt,order=2),
                probe.acf("y",lags=c(0,3),transform=sqrt),
                mean=probe.mean("y",transform=sqrt)
                ),
              transform=TRUE,
              nsim=1000,
              seed=1066L
              )
plot(pb)

## ----ricker-probe-match-calc,eval=F,results="markup"---------------------
## pm <- probe.match(
##                   pb.guess,
##                   est=c("r","sigma","phi"),
##                   transform=TRUE,
##                   method="Nelder-Mead",
##                   maxit=2000,
##                   seed=1066L,
##                   reltol=1e-8
##                   )

## ----ricker-probe.match-eval,echo=F,eval=T,results="hide",cache=F--------
bake("ricker-probe-match.rda",{
  pm <- probe.match(
                    pb.guess,
                    est=c("r","sigma","phi"),
                    transform=TRUE,
                    method="Nelder-Mead",
                    maxit=2000,
                    seed=1066L,
                    reltol=1e-8
                    )
})

## ----ricker-mif-calc,eval=FALSE,tidy=FALSE-------------------------------
## mf <- mif(ricker, start = guess, Nmif = 100, Np = 1000, transform = TRUE,
##           cooling.fraction = 0.95^50, var.factor = 2, ic.lag = 3,
##           rw.sd=c(r = 0.1, sigma = 0.1, phi = 0.1), max.fail = 50)
## mf <- continue(mf, Nmif = 500, max.fail = 20)

## ----ricker-mif-eval,echo=F,eval=T,cache=F,results="hide"----------------
bake("ricker-mif.rda",{
  save.seed <- .Random.seed
  set.seed(718086921L)
  mf <- mif(ricker, start = guess, Nmif = 100, Np = 1000, transform = TRUE,
            cooling.fraction = 0.95^50, var.factor = 2, ic.lag = 3,
            rw.sd=c(r = 0.1, sigma = 0.1, phi = 0.1), max.fail = 50)
  mf <- continue(mf, Nmif = 500, max.fail = 20)
  .Random.seed <<- save.seed
})

## ----ricker-comparison,eval=T,echo=F,cache=F-----------------------------
bake("ricker-comparison.rda",{
  require(plyr)
  require(magrittr)
  require(foreach)
  require(doMC)
  registerDoMC()

  rbind(guess=guess,
        truth=coef(ricker),
        MLE=coef(mf),
        MSLE=coef(pm)) -> comp

  save.seed <- .Random.seed
  set.seed(1182206495L,kind="L'Ecuyer")

  foreach (theta=iter(comp,"row"),.combine=rbind) %do% {

    foreach (i=seq_len(10),.combine=c,.inorder=FALSE,
             .options.multicore=list(set.seed=TRUE)) %dopar% {
               logLik(pfilter(ricker,params=theta[1,],Np=10000))
             } %>% logmeanexp(se=TRUE) -> pf

    foreach (i=seq_len(10),.combine=c,.inorder=FALSE,
             .options.multicore=list(set.seed=TRUE)) %dopar% {
               logLik(probe(ricker,params=theta[1,],nsim=10000,
                            probes=plist))
             } %>% logmeanexp(se=TRUE) -> pb

    cbind(theta,loglik=pf[1],loglik.se=pf[2],logsynlik=pb[1],logsynlik.se=pb[2])
  } %>% subset(select=-c(N.0,e.0)) -> comp

  .Random.seed <<- save.seed
})

## ----ricker-comparison-show,echo=F,results="asis"------------------------
require(xtable)
colnames(comp) <- c("$r$","$\\sigma$","$\\phi$",
                    "$\\loglikMC$","s.e.($\\loglikMC$)",
                    "$\\synloglikMC$","s.e.($\\synloglikMC$)")
print(xtable(comp,align="r|ccccccc",digits=c(0,1,3,1,1,2,1,2)))

## ----abc-load,eval=F,echo=F,tidy=F---------------------------------------
## plist <- list(probe.mean(var = "Y", transform = sqrt),
##               probe.acf("Y", lags = c(0, 5, 10, 20)),
##               probe.marginal("Y", ref = obs(gompertz)))
## psim <- probe(gompertz, probes = plist, nsim = 500)
## scale.dat <- apply(psim$simvals, 2, sd)

## ----abc-inner,eval=FALSE,echo=FALSE-------------------------------------
## abc(pomp(gompertz, dprior = gompertz.dprior), Nabc = 4e6,
##     probes = plist, epsilon = 2, scale = scale.dat,
##     proposal=mvn.diag.rw(c(r = 0.01, sigma = 0.01, tau = 0.01)))

## ----abc-demo,eval=FALSE,tidy=FALSE--------------------------------------
## plist <- list(probe.mean(var = "Y", transform = sqrt),
##               probe.acf("Y", lags = c(0, 5, 10, 20)),
##               probe.marginal("Y", ref = obs(gompertz)))
## psim <- probe(gompertz, probes = plist, nsim = 500)
## scale.dat <- apply(psim$simvals, 2, sd)
## abc1 <- foreach(i = 1:5, .combine = c) %dopar% {
##   abc(pomp(gompertz, dprior = gompertz.dprior), Nabc = 4e6,
##       probes = plist, epsilon = 2, scale = scale.dat,
##       proposal=mvn.diag.rw(c(r = 0.01, sigma = 0.01, tau = 0.01)))
## }

## ----abc-eval,echo=F,results="hide",cache=F------------------------------
require(pomp)

plist <- list(probe.mean(var = "Y", transform = sqrt),
              probe.acf("Y", lags = c(0, 5, 10, 20)),
              probe.marginal("Y", ref = obs(gompertz)))
psim <- probe(gompertz, probes = plist, nsim = 500)
scale.dat <- apply(psim$simvals, 2, sd)

bake("abc.rda",{

  pompExample(gompertz)

  tic <- Sys.time()
  require(doMC)
  require(foreach)
  registerDoMC()

  save.seed <- .Random.seed
  set.seed(334388458L,kind="L'Ecuyer")

  abc1 <- foreach(
    i=1:5,
    .inorder=FALSE,
    .packages="pomp",
    .combine=c,
    .options.multicore=list(set.seed=TRUE)
  ) %dopar% {
    abc(pomp(gompertz, dprior = gompertz.dprior), Nabc = 4e6,
        probes = plist, epsilon = 2, scale = scale.dat,
        proposal=mvn.diag.rw(c(r = 0.01, sigma = 0.01, tau = 0.01)))
  }

toc <- Sys.time()
abcTime <- toc-tic
.Random.seed <<- save.seed

abc.traces <- conv.rec(abc1,c("r","sigma","tau"))
abc.traces <- window(abc.traces,start=2000001,thin=400)
ess.abc <- effectiveSize(abc.traces)
rm(abc1,save.seed,tic,toc)
})


## ----abc-diagnostics,results="hide",fig.show="hide",echo=FALSE,eval=TRUE----
gelman.diag(abc.traces)
gelman.plot(abc.traces)
autocorr.plot(abc.traces[[1]])
hist(rle(unlist(abc.traces[,"r"]))$length)

## ----abc-pmmc-compare,echo=F,fig.width=7,fig.height=3,cache=T------------
require(ggplot2)
require(grid)
require(plyr)
require(reshape2)
require(magrittr)

ldply(list(pmcmc=ldply(pmcmc.traces),abc=ldply(abc.traces)),.id='method') %>%
  melt(id="method") %>%
  mutate(log.value=log10(value)) -> traces

coef(gompertz,c("r","sigma","tau")) %>% as.list() %>% as.data.frame() %>%
  melt(id=NULL) %>%
  mutate(log.value=log10(value)) -> truth

traces %>%
  ggplot(mapping=aes(x=log.value,linetype=method))+
  geom_density(adjust=5)+
  geom_vline(data=truth,mapping=aes(xintercept=log.value))+
  scale_linetype_manual(values=c(abc=2,pmcmc=1))+
  facet_grid(~variable,scales="free_x",
             labeller=labeller(variable=list(sigma=expression(log[10]~sigma),
                                             tau=expression(log[10]~tau),
                                             r=expression(log[10]~r))))+
  labs(x="",y="",linetype="")+
  theme_classic()+
  theme(legend.position=c(0.35,0.7),
        strip.background=element_rect(fill=NA,color=NA),
        strip.text=element_text(size=12),
        panel.margin=unit(4,"mm"))

## ----first-nlf,eval=F,tidy=F---------------------------------------------
## nlf1 <- nlf(gompertz, nasymp = 1000, nconverge = 1000, lags = c(2, 3),
##             start = c(r = 1, K = 2, sigma = 0.5, tau = 0.5, X.0 = 1),
##             est = c("r", "sigma", "tau"), transform = TRUE)

## ----nlf-mif-comp-setup,eval=T,echo=F,results="hide"---------------------
pompExample(gompertz)
set.seed(4897341L)
## number of replicates:
R <- 10
pompExample(gompertz)
estpars <- c("r","sigma","tau")
gompList <- simulate(gompertz,nsim=R)


## ----nlf-mif-compare-eval,echo=F,eval=T,results="hide"-------------------
bake("nlf-mif-compare.rda",{
  require(doMC)
  require(foreach)
  registerDoMC()

  save.seed <- .Random.seed
  set.seed(816326853L,kind="L'Ecuyer")

  tic <- Sys.time()
  cmp1 <- foreach(gomp=gompList,
                  .inorder=FALSE,.packages="pomp",.combine=rbind,
                  .options.multicore=list(set.seed=TRUE)
  ) %dopar% {
    true.lik <- pfilter(gomp,Np=10000)
    true.sql <- nlf(
      gomp,
      nasymp=10000,
      eval.only=TRUE,
      lags=c(2,3)
    )

  ## start at the truth:
  theta.guess <- coef(gomp)

  tic <- Sys.time()
  mif1 <- mif(gomp,Nmif=100,start=theta.guess,transform=TRUE,
    rw.sd=c(r=0.02,sigma=0.02,tau=0.05),Np=1000,var.factor=4,ic.lag=10,
    cooling.type="geometric",cooling.fraction=0.5)
  mif.lik <- pfilter(mif1,Np=10000)
  toc <- Sys.time()
  mif.time <- toc-tic
  units(mif.time) <- "secs"
  mif.sql <- nlf(mif1,nasymp=10000,eval.only=TRUE,lags=c(2,3))

  tic <- Sys.time()
  nlf1 <- nlf(gomp,start=theta.guess,transform=TRUE,est=estpars,lags=c(2,3))
  toc <- Sys.time()
  nlf.time <- toc-tic
  units(nlf.time) <- "secs"
  nlf.lik <- pfilter(nlf1,Np=10000)
  nlf.sql <- nlf(nlf1,nasymp=10000,eval.only=TRUE)

  c(
    trueLik=logLik(true.lik),
    trueSQL=logLik(true.sql),
    mifLik=logLik(mif.lik),
    mifSQL=logLik(mif.sql),
    nlfLik=logLik(nlf.lik),
    nlfSQL=logLik(nlf.sql),
    mifTime=mif.time,
    nlfTime=nlf.time
  )
}
toc <- Sys.time()
nlf.mif.time <- toc-tic
.Random.seed <<- save.seed
cmp1 <- as.data.frame(cmp1)
})

## ----nlf-mif-plot,echo=F,fig.width=8,fig.height=3.5----------------------
require(ggplot2)
require(grid)
require(plyr)
require(reshape2)
require(magrittr)

plA <- cmp1 %>%
  ggplot(mapping=aes(y=mifLik-trueLik,x=nlfLik-trueLik))+
  geom_point()+
  geom_abline(slope=1,yintercept=0,linetype=3)+
  geom_hline(yintercept=0,linetype=3)+
  geom_vline(xintercept=0,linetype=3)+
  expand_limits(x=c(-1,1),y=c(-1,1))+
  labs(y=expression(hat("\u2113")(hat(theta))-hat("\u2113")(theta)),
       x=expression(hat("\u2113")(tilde(theta))-hat("\u2113")(theta)))+
  theme_classic()

plB <- cmp1 %>%
  ggplot(mapping=aes(y=mifSQL-trueSQL,x=nlfSQL-trueSQL))+
  geom_point()+
  geom_abline(slope=1,yintercept=0,linetype=3)+
  geom_hline(yintercept=0,linetype=3)+
  geom_vline(xintercept=0,linetype=3)+
  expand_limits(x=c(-1,1),y=c(-1,1))+
  labs(y=expression(hat("\u2113")[Q](hat(theta))-hat("\u2113")[Q](theta)),
       x=expression(hat("\u2113")[Q](tilde(theta))-hat("\u2113")[Q](theta)))+
  theme_classic()

grid.newpage()
pushViewport(viewport(layout=grid.layout(1,2)))
print(plA,vp=viewport(layout.pos.row=1,layout.pos.col=1))
print(plB,vp=viewport(layout.pos.row=1,layout.pos.col=2))
grid.text("A",x=unit(0.1,"npc"),y=unit(1,"npc"),
          vp=viewport(layout.pos.row=1,layout.pos.col=1),
          hjust=0.5,vjust=1,
          gp=gpar(fontsize=16,fontface="bold"))
grid.text("B",x=unit(0.1,"npc"),y=unit(1,"npc"),
          vp=viewport(layout.pos.row=1,layout.pos.col=2),
          hjust=0.5,vjust=1,
          gp=gpar(fontsize=16,fontface="bold"))
popViewport()

## ----sir-measmodel,tidy=F------------------------------------------------
rmeas <- "
  cases = rnbinom_mu(theta, rho * H);
"
dmeas <- "
  lik = dnbinom_mu(cases, theta, rho * H, give_log);
"

## ----sir-proc-sim-def,tidy=F---------------------------------------------
sir.step <- "
  double rate[6];
  double dN[6];
  double P;
  P = S + I + R;
  rate[0] = mu * P;       // birth
  rate[1] = beta * I / P; // transmission
  rate[2] = mu;           // death from S
  rate[3] = gamma;        // recovery
  rate[4] = mu;           // death from I
  rate[5] = mu;           // death from R
  dN[0] = rpois(rate[0] * dt);
  reulermultinom(2, S, &rate[1], dt, &dN[1]);
  reulermultinom(2, I, &rate[3], dt, &dN[3]);
  reulermultinom(1, R, &rate[5], dt, &dN[5]);
  S += dN[0] - dN[1] - dN[2];
  I += dN[1] - dN[3] - dN[4];
  R += dN[3] - dN[5];
  H += dN[1];
"

## ----sir-pomp-def,eval=T,echo=T,results="hide",tidy=F--------------------
pomp(data = data.frame(cases = NA, time = seq(0, 10, by=1/52)),
     times = "time", t0 = -1/52, dmeasure = Csnippet(dmeas),
     rmeasure = Csnippet(rmeas), rprocess = euler.sim(
       step.fun = Csnippet(sir.step), delta.t = 1/52/20),
     statenames = c("S", "I", "R", "H"),
     paramnames = c("gamma", "mu", "theta", "beta", "popsize",
       "rho", "S.0", "I.0", "R.0"), zeronames=c("H"),
     initializer=function(params, t0, ...) {
       fracs <- params[c("S.0", "I.0", "R.0")]
       setNames(c(round(params["popsize"]*fracs/sum(fracs)),0),
                c("S","I","R","H"))
     }, params = c(popsize = 500000, beta = 400, gamma = 26,
          mu = 1/50, rho = 0.1, theta = 100, S.0 = 26/400,
          I.0 = 0.002, R.0 = 1)) -> sir1
simulate(sir1, seed = 1914679908L) -> sir1

## ----sir1-plot,echo=F,fig.height=5---------------------------------------
ops <- options(scipen=-10)
plot(sir1,mar=c(0,5,2,0))
options(ops)

## ----birthdat,eval=T,echo=F----------------------------------------------
birthdat <- data.frame(time=seq(-1,11,by=1/12))
birthdat$births <- 5e5*bspline.basis(birthdat$time,nbasis=5)%*%c(0.018,0.019,0.021,0.019,0.015)
save.seed <- .Random.seed
set.seed(5853712L)
birthdat$births <- ceiling(rlnorm(
                                  n=nrow(birthdat),
                                  meanlog=log(birthdat$births),
                                  sdlog=0.001
                                  ))
.Random.seed <<- save.seed

## ----complex-sir-def,echo=T,eval=T,results="hide",tidy=F-----------------
seas.sir.step <- "
  double rate[6];
  double dN[6];
  double Beta;
  double dW;
  Beta = exp(b1 + b2 * cos(M_2PI * Phi) + b3 * sin(M_2PI * Phi));
  rate[0] = births;                // birth
  rate[1] = Beta * (I + iota) / P; // infection
  rate[2] = mu;                    // death from S
  rate[3] = gamma;                 // recovery
  rate[4] = mu;                    // death from I
  rate[5] = mu;                    // death from R
  dN[0] = rpois(rate[0] * dt);
  reulermultinom(2, S, &rate[1], dt, &dN[1]);
  reulermultinom(2, I, &rate[3], dt, &dN[3]);
  reulermultinom(1, R, &rate[5], dt, &dN[5]);
  dW = rnorm(dt, sigma * sqrt(dt));
  S += dN[0] - dN[1] - dN[2];
  I += dN[1] - dN[3] - dN[4];
  R += dN[3] - dN[5];
  P = S + I + R;
  Phi += dW;
  H += dN[1];
  noise += (dW - dt) / sigma;
"
pomp(sir1, rprocess = euler.sim(
  step.fun = Csnippet(seas.sir.step), delta.t = 1/52/20),
  dmeasure = Csnippet(dmeas), rmeasure = Csnippet(rmeas),
  covar = birthdat, tcovar = "time", zeronames = c("H", "noise"),
  statenames = c("S", "I", "R", "H", "P", "Phi", "noise"),
  paramnames = c("gamma", "mu", "popsize", "rho","theta","sigma",
                 "S.0", "I.0", "R.0", "b1", "b2", "b3", "iota"),
  initializer = function(params, t0, ...) {
    fracs <- params[c("S.0", "I.0", "R.0")]
    setNames(c(round(params["popsize"]*c(fracs/sum(fracs),1)),0,0,0),
             c("S","I","R","P","H","Phi","noise"))
  }, params = c(popsize = 500000, iota = 5, b1 = 6, b2 = 0.2,
                b3 = -0.1, gamma = 26, mu = 1/50, rho = 0.1, theta = 100,
                sigma = 0.3, S.0 = 0.055, I.0 = 0.002, R.0 = 0.94)) -> sir2
simulate(sir2, seed = 619552910L) -> sir2

## ----sir2-plot,echo=F,fig.height=6.5-------------------------------------
ops <- options(scipen=-10)
plot(sir2,mar=c(0,5,2,0))
options(ops)

## ----timing2,cache=F-----------------------------------------------------
bigtock <- Sys.time()
totalSweaveTime <- bigtock-bigtick

if (file.exists("timing.rda")) {
  load("timing.rda")
} else {
  save(totalSweaveTime,file='timing.rda',compress='xz')
}

print(Sys.info())
print(sessionInfo())

print(mifTime)
print(pmcmcTime)
print(abcTime)
print(nlf.mif.time)
print(totalSweaveTime)

