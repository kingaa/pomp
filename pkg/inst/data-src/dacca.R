require(pomp)

cholera <- dget("dacca-cholera.q")
census <- dget("dacca-census.q")
mle <- dget("dacca-mle.q")

covar.dt <- 0.01
nbasis <- as.integer(mle["nbasis"])
nrstage <- as.integer(mle["nrstage"])

t0 <- with(cholera,2*time[1]-time[2])
tcovar <- seq(from=t0,to=max(cholera$time)+2/12,by=covar.dt)
covartable <- data.frame(
                         time=tcovar,
                         seas=periodic.bspline.basis(tcovar-1/12,nbasis=nbasis,degree=3,period=1),
                         pop=predict(smooth.spline(x=census$year,y=census$census),x=tcovar)$y,
                         dpopdt=predict(smooth.spline(x=census$year,y=census$census),x=tcovar,deriv=1)$y,
                         trend=tcovar-mean(tcovar)
                         )

pomp(
     data=cholera[c("cholera.deaths","time")],
     times='time',
     t0=t0,
     delta.t=1/240,
     nrstage = nrstage,
     rprocess = euler.simulate,
     rmeasure = "_cholmodel_norm_rmeasure",
     dmeasure = "_cholmodel_norm_dmeasure",
     step.fun = "_cholmodel_one",
     covar=covartable,
     tcovar='time',
     statenames = c("S","I","Rs","R1","M","W","count"),
     paramnames = c("tau","gamma","eps","delta","deltaI",
       "omega1","sd.beta","beta.trend","log.beta1",
       "alpha","rho","clin","nbasis","nrstage"),
     covarnames = c("pop","dpopdt","seas.1","trend"),
     zeronames = c("M","count"),
     PACKAGE = "pomp",
     initializer = function (params, t0, covars, all.state.names, comp.names, nrstage, ...) {
       all.state.names <- c("S","I","Rs","R1","R2","R3","M","W","count")
       comp.names <- c("S","I","Rs",paste("R",1:nrstage,sep=''))
       comp.ic.names <- paste(comp.names,"0",sep='.')
       states <- numeric(length(all.state.names))
       names(states) <- all.state.names
       ## translate fractions into initial conditions
       frac <- exp(params[comp.ic.names])
       states[comp.names] <- round(covars['pop']*frac/sum(frac))
       states
     }
     ) -> dacca

## parameter transformations for fitting
cholera.transform <- function (params, dir = c("forward","inverse")) {
  dir <- match.arg(dir)
  r <- length(dim(params))
  r <- length(dim(params))
  log.trans=c(                 # parameters to log transform
    "gamma","eps","rho","delta","deltaI","alpha",
    "tau","sd.beta",
    paste("omega",seq(length=nbasis),sep=''),
    "S.0","I.0","Rs.0","R1.0","R2.0","R3.0"
    )
  logit.trans="clin"
  logit <- function(p){log(p/(1-p))}      # (0,1) -> (-inf,inf)
  expit <- function(r){1/(1+exp(-r))}     # (-inf,inf) -> (0,1)
  par.trans.vec <- function (x) {
    x[logit.trans] <- logit(x[logit.trans])
    x[log.trans] <- log(x[log.trans])
    x
  }
  par.untrans.vec <- function (x) {
    x[logit.trans] <- expit(x[logit.trans])
    x[log.trans] <- exp(x[log.trans])
    x
  }
  switch(
         dir,
         forward=if (r > 1) apply(params,2:r,par.trans.vec) else par.trans.vec(params),
         inverse=if (r > 1) apply(params,2:r,par.untrans.vec) else par.untrans.vec(params)
         )
}

coef(dacca) <- cholera.transform(mle,dir="forward")
