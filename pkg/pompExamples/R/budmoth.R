budmoth.sim <- function (which) {
  if (missing(which)) {
    cat("available datasets:",
        sQuote(c("food","para1","para2","tri")),"\n")
  } else {
    w <- as.character(substitute(which))
    simulate(
             pomp(
                  data=data.frame(
                    time=seq(from=0,to=60,by=1),
                    Qobs=NA,Nobs=NA,Sobs=NA
                    ),
                  time="time",
                  t0=-1,
                  params=switch(
                    w,
                    tri=c(
                      alpha=0.5, sig.alpha=0.1, gam=50, lambda=22, sig.lambda=0.25, g=0.08, delta=10,
                      a=1.7, sig.a=0.1, w=0.15, beta0=0, beta1=35, u=0.9,
                      sigQobs=0.03, sigNobs=0.5, sigSobs=0.1,
                      Q.0=0.96, N.0=0.02, S.0=0.22
                      ),
                    food=c(
                      alpha=0.5, sig.alpha=0.1, gam=20, lambda=5, sig.lambda=0.25, g=0.02, delta=10,
                      a=1, sig.a=0.1, w=0, beta0=0, beta1=35, u=0.9,
                      sigQobs=0.03, sigNobs=0.5, sigSobs=0.1,
                      Q.0=0.96, N.0=0.02, S.0=0.22
                      ),
                    para1=c(
                      alpha=0.5, sig.alpha=0.1, gam=50, lambda=22, sig.lambda=0.25, g=0.08, delta=0.5,
                      a=1.7, sig.a=0.1, w=0.15, beta0=0, beta1=35, u=0.9,
                      sigQobs=0.03, sigNobs=0.5, sigSobs=0.1,
                      Q.0=0.96, N.0=0.02, S.0=0.22
                      ),
                    para2=c(
                      alpha=0.5, sig.alpha=0.1, gam=50, lambda=10, sig.lambda=5, g=0.08, delta=0.5,
                      a=1.7, sig.a=1, w=0.15, beta0=0, beta1=35, u=0.9,
                      sigQobs=0.03, sigNobs=0.5, sigSobs=0.1,
                      Q.0=0.96, N.0=0.02, S.0=0.22
                      ),
                    stop("unrecognized dataset ",sQuote(w),call.=FALSE)
                    ),
                  rprocess=euler.sim(
                    step.fun="budmoth_map",
                    delta.t=1,
                    PACKAGE="pompExamples"
                    ),
                  dprocess=onestep.dens(
                    dens.fun="budmoth_density",
                    PACKAGE="pompExamples"
                    ),
                  rmeasure="budmoth_rmeasure",
                  dmeasure="budmoth_dmeasure",
                  skeleton.type="map",
                  skeleton="budmoth_skeleton",
                  PACKAGE="pompExamples",
                  paramnames=c(
                    "alpha","sig.alpha","gam","lambda","sig.lambda",
                    "g","delta","a","sig.a",
                    "w","beta0","beta1","u",
                    "sigQobs","sigNobs","sigSobs"
                    ),
                  statenames=c(
                    "Alpha","Lambda","A","Q","N","S"
                    ),
                  obsnames=c("Qobs","Nobs","Sobs"),
                  initializer=function (params, t0, ...) {
                    x <- c(params[c("Q.0","N.0","S.0")],c(0,0,0))
                    names(x) <- c("Q","N","S","Alpha","Lambda","A")
                    x
                  },
                  logitvar=c("alpha","Q.0","S.0","u"),
                  logvar=c(
                    "sig.alpha","gam","lambda","sig.lambda",
                    "g","delta","a","w","sig.a","beta1","sigQobs",
                    "sigNobs", "sigSobs","N.0"
                    ),
                  parameter.transform=function (params, logitvar, logvar, ...) {
                    params[logitvar] <- plogis(params[logitvar])
                    params[logvar] <- exp(params[logvar])
                    params
                  },
                  parameter.inv.transform=function (params, logitvar, logvar, ...) {
                    params[logitvar] <- qlogis(params[logitvar])
                    params[logvar] <- log(params[logvar])
                    params
                  }
                  ),
             seed=switch(
               w,
               tri=1691699385L,
               food=1054866677L,
               para1=1116757478L,
               para2=1361101458L,
               stop("unrecognized dataset ",sQuote(w),call.=FALSE)
               )
             )
  }
}
