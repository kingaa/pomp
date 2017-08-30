## ----opts,include=FALSE,cache=FALSE--------------------------------------
options(stringsAsFactors=FALSE)
library(ggplot2)
theme_set(theme_bw())
set.seed(2028866059L)

## ----parus-data----------------------------------------------------------
loc <- url("http://kingaa.github.io/pomp/vignettes/parus.csv")
dat <- read.csv(loc)
head(dat)
plot(pop~year,data=dat,type='o')

## ----parus-pomp1---------------------------------------------------------
library(pomp)
parus <- pomp(dat,times="year",t0=1959)

## ----parus-plot1---------------------------------------------------------
plot(parus)

## ----parus-sim-defn------------------------------------------------------
stochStep <- Csnippet("
  N = r*N*exp(-c*N+rnorm(0,sigma));
")

pomp(
  parus,
  rprocess=discrete.time.sim(step.fun=stochStep,delta.t=1),
  initializer=Csnippet("N = N_0;"),
  paramnames=c("r","c","sigma","N_0"),
  statenames=c("N")
) -> parus

## ----ricker-first-sim----------------------------------------------------
sim <- simulate(parus, params=c(N_0=1,r=12,c=1,sigma=0.5),
  as.data.frame=TRUE, states=TRUE)

plot(N~time,data=sim,type='o')

## ----parus-rmeas-defn----------------------------------------------------
rmeas <- Csnippet("pop = rpois(phi*N);")

## ----parus-dmeas-defn----------------------------------------------------
dmeas <- Csnippet("lik = dpois(pop,phi*N,give_log);")

## ----parus-add-meas------------------------------------------------------
pomp(parus,
     rmeasure=rmeas,
     dmeasure=dmeas,
     statenames=c("N"),
     paramnames=c("phi")
) -> parus

## ----ricker-add-params---------------------------------------------------
coef(parus) <- c(N_0=1,r=20,c=1,sigma=0.1,phi=200)

## ----ricker-second-sim,results='markup',fig.height=6,fig.width=5---------
library(ggplot2)

sims <- simulate(parus,nsim=3,as.data.frame=TRUE,include.data=TRUE)

ggplot(data=sims,
       mapping=aes(x=time,y=pop))+
  geom_line()+
  facet_wrap(~sim,ncol=1,scales="free_y")

## ----plot-ricker---------------------------------------------------------
plot(parus)

## ----sim-ricker1---------------------------------------------------------
x <- simulate(parus)

## ------------------------------------------------------------------------
class(x)
plot(x)

## ------------------------------------------------------------------------
y <- as.data.frame(parus)
head(y)
head(simulate(parus,as.data.frame=TRUE))

## ------------------------------------------------------------------------
x <- simulate(parus,nsim=10)
class(x)
sapply(x,class)
x <- simulate(parus,nsim=10,as.data.frame=TRUE)
head(x)
str(x)

## ----fig.height=8--------------------------------------------------------
library(ggplot2)
x <- simulate(parus,nsim=9,as.data.frame=TRUE,include.data=TRUE)
ggplot(data=x,aes(x=time,y=pop,group=sim,color=(sim=="data")))+
  geom_line()+guides(color=FALSE)+
  facet_wrap(~sim,ncol=2)

## ----coef-ricker---------------------------------------------------------
coef(parus)

## ------------------------------------------------------------------------
theta <- coef(parus)
theta[c("r","N_0")] <- c(5,3)

x <- simulate(parus,params=theta)

plot(x,var="pop")

## ------------------------------------------------------------------------
coef(parus,c("r","N_0","sigma")) <- c(5,1.5,0.1)
coef(parus)
plot(simulate(parus),var=c("pop","N"))

## ----ricker-pfilter------------------------------------------------------
pf <- pfilter(parus,Np=1000)
class(pf)
logLik(pf)

plot(pf,var=c("ess","cond.loglik"))

## ----sir-measmodel-------------------------------------------------------
rmeas <- "
  cases = rnbinom_mu(theta, rho * H);
"

dmeas <- "
  lik = dnbinom_mu(cases, theta, rho * H, give_log);
"

## ----sir-proc-sim-def----------------------------------------------------
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

## ----sir-pomp-def,eval=T,echo=T,results="hide"---------------------------
sir1 <- pomp(
  data = data.frame(cases = NA, time = seq(0, 10, by = 1/52)),
  times = "time", t0 = -1/52, 
  dmeasure = Csnippet(dmeas),
  rmeasure = Csnippet(rmeas), 
  rprocess = euler.sim(step.fun = Csnippet(sir.step), delta.t = 1/52/20),
  statenames = c("S", "I", "R", "H"),
  paramnames = c("gamma", "mu", "theta", "beta", "popsize",
                 "rho", "S.0", "I.0", "R.0"), 
  zeronames = "H",
  initializer = Csnippet("
    double sum = S_0 + I_0 + R_0;
    S = nearbyint(popsize * S_0 / sum);
    I = nearbyint(popsize * I_0 / sum);
    R = nearbyint(popsize * R_0 / sum);
    H = 0;
    "),
  params = c(popsize = 500000, beta = 400, gamma = 26,
             mu = 1/50, rho = 0.1, theta = 100, S.0 = 26/400,
             I.0 = 0.002, R.0 = 1))

sir1 <- simulate(sir1, seed = 1914679908L)

## ----fig.height=5,echo=FALSE---------------------------------------------
ops <- options(scipen=-10)
plot(sir1,mar=c(0,5,2,0))
options(ops)

## ----birthdat,eval=T,echo=F,results="hide"-------------------------------
##' Construct some fake birthrate data.
birthdat <- data.frame(time=seq(-1,11,by=1/12))
birthdat$births <- 5e5*bspline.basis(birthdat$time,nbasis=5)%*%c(0.018,0.019,0.021,0.019,0.015)
freeze(seed=5853712L,{
  birthdat$births <- ceiling(rlnorm(
                                    n=nrow(birthdat),
                                    meanlog=log(birthdat$births),
                                    sdlog=0.001
                                    ))
})

## ----complex-sir-def,echo=T,eval=T,results="hide"------------------------
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

sir2 <- pomp(
  sir1, 
  rprocess = euler.sim(
    step.fun = Csnippet(seas.sir.step), delta.t = 1/52/20
  ),
  dmeasure = Csnippet(dmeas), 
  rmeasure = Csnippet(rmeas),
  covar = birthdat, tcovar = "time", 
  zeronames = c("H", "noise"),
  statenames = c("S", "I", "R", "H", "P", "Phi", "noise"),
  paramnames = c("gamma", "mu", "popsize", "rho", "theta", 
                 "sigma", "S.0", "I.0", "R.0", 
                 "b1", "b2", "b3", "iota"),
  initializer = Csnippet("
    double sum = S_0 + I_0 + R_0;
    S = nearbyint(popsize * S_0 / sum);
    I = nearbyint(popsize * I_0 / sum);
    R = nearbyint(popsize * R_0 / sum);
    P = S + I + R;
    H = 0;
    Phi = 0;
    noise = 0;
    "),
  params = c(popsize = 500000, iota = 5,
             b1 = 6, b2 = 0.2, b3 = -0.1,
             gamma = 26, mu = 1/50, rho = 0.1, theta = 100,
             sigma = 0.3, S.0 = 0.055, I.0 = 0.002, R.0 = 0.94)
  )

sir2 <- simulate(sir2, seed = 619552910L)

## ----sir2-plot,echo=F,fig.height=6.5-------------------------------------
ops <- options(scipen=-10)
plot(sir2,mar=c(0,5,2,0))
options(ops)

## ----parus-pfilter1,warning=T--------------------------------------------
pf <- replicate(10, pfilter(parus,Np=5000))
plot(pf[[1]])
ll <- sapply(pf,logLik)
logmeanexp(ll,se=TRUE)

## ----parus-mif1----------------------------------------------------------
mif2(parus, Nmif=30, Np=1000, 
     cooling.fraction.50=0.8,cooling.type="geometric",
     rw.sd=rw.sd(r=0.02,sigma=0.02,phi=0.02,N_0=ivp(0.1))
) -> mf

plot(mf)

## ----parus-pfilter2------------------------------------------------------
pf <- replicate(5, pfilter(mf,Np=1000))
ll <- sapply(pf,logLik)
logmeanexp(ll,se=TRUE)

## ----parus-sim3----------------------------------------------------------
simulate(mf,nsim=5,as.data.frame=TRUE,include.data=TRUE) -> sims

library(ggplot2)
library(plyr)
ggplot(sims, 
       mapping=aes(x=time,y=pop,group=sim,color=sim=="data"))+
  geom_line()+
  guides(color=FALSE)+
  facet_wrap(~sim)

## ----parus-probe1,fig.height=6,fig.width=6-------------------------------
probe(mf, nsim=200, 
      probes=list(mean=probe.mean("pop"),
                  sd=probe.sd("pop"),
                  probe.acf("pop",transform=sqrt,lags=c(1,2)),
                  probe.quantile("pop",prob=c(0.2,0.8))
      )) -> pb
      
plot(pb)

## ----parus-probematch1---------------------------------------------------
probe.match(pb, nsim=200, method="subplex",
            seed = 669237763L,
            est = c("N_0","r")) -> pm
summary(pm)

## ----parus-pmcmc1--------------------------------------------------------
priorDens <- "
  lik = dnorm(sigma,0.2,1,1)+
    dnorm(phi,200,100,1)+
    dexp(r,0.1,1);
  if (!give_log) lik = exp(lik);
"

pmcmc(pomp(mf, dprior=Csnippet(priorDens),
           paramnames=c("sigma","phi","r")),
      Nmcmc = 500, Np = 1000, 
      proposal = mvn.rw.adaptive(
        rw.sd=c(N_0=0.1, sigma=0.02, r=0.02, phi=0.02)
      )) -> pmh

plot(pmh)

## ----parus-coda1---------------------------------------------------------
library(coda)

conv.rec(pmh) -> trace
class(trace)

autocorr.diag(trace[,c("r","sigma","phi")])

## ------------------------------------------------------------------------
pfilter <-  function (object, params, Np,
                      tol = 1e-17, warn = TRUE,
                      max.fail = 0,
                      ...) {
  if (missing(Np))
    Np <- NCOL(params)
  ntimes <- length(time(object))
  if (is.null(dim(params))) {
    params <- matrix(params,nrow=length(params),ncol=Np,
                     dimnames=list(names(params),NULL))
  }
  npars <- nrow(params)
  paramnames <- rownames(params)
  if (is.null(paramnames))
    stop("pfilter error: 'params' must have rownames")
  
  ## sample from the X0 distribution
  xstart <- init.state(object,params=params,t0=object@t0)
  statenames <- rownames(xstart)
  nvars <- nrow(xstart)
  
  loglik <- rep(NA,ntimes)
  eff.sample.size <- rep(NA,ntimes)
  nfail <- 0
  
  times <- c(object@t0,time(object))
  x <- xstart
  
  for (nt in seq(length=ntimes)) {
    
    ## advance the state variables according to the process model
    X <- try(
      rprocess(object,x=x,times=times[c(nt,nt+1)],
               params=params)[,,2,drop=FALSE],
      silent=FALSE
    )
    if (inherits(X,'try-error'))
      stop("pfilter error: process simulation error")
    
    x[,] <- X                 # drop the third dimension
    
    ## determine the weights
    weights <- try(
      dmeasure(object,y=object@data[,nt,drop=FALSE],x=X,
               times=times[nt+1],params=params),
      silent=FALSE
    )
    if (inherits(weights,'try-error'))
      stop("pfilter error: error in calculation of weights")
    if (any(is.na(weights))) {
      stop("pfilter error: dmeasure returns NA")
    }
    
    ## test for failure to filter
    dim(weights) <- NULL
    failures <- weights < tol
    all.fail <- all(failures)
    if (all.fail) {                     # all particles are lost
      if (warn)
        message("filtering failure at time t = ",times[nt+1])
      nfail <- nfail+1
      if (nfail > max.fail)
        stop('pfilter error: too many filtering failures')
      loglik[nt] <- log(tol)          # worst log-likelihood
      weights <- rep(1/Np,Np)
      eff.sample.size[nt] <- 0
    } else {                  # not all particles are lost
      ## compute log-likelihood
      loglik[nt] <- log(mean(weights))  
      weights[failures] <- 0
      weights <- weights/sum(weights)
      ## compute effective sample-size
      eff.sample.size[nt] <- 1/(weights%*%weights) 
    }
    
     if (!all.fail) {
      sample <- .Call(systematic_resampling,weights)
      x <- x[,sample,drop=FALSE]
      params <- params[,sample,drop=FALSE]
    }
  }
  
  list(
    pred.mean=pred.m,
    pred.var=pred.v,
    filter.mean=filt.m,
    eff.sample.size=eff.sample.size,
    cond.loglik=loglik,
    nfail=nfail,
    loglik=sum(loglik)
  )
}

