require(pomp.devel)

params <- read.csv("budmoth-params.csv",sep=";")
datasets <- rownames(params)
params <- params[!(names(params)%in%c("seed"))]

load("budmoth-simdata.rda")

budmoth.sim <- list()
for (d in datasets) {
  po <- pomp(
             data=subset(budmoth.sim.data,dataset==d,select=c(time,Qobs,Nobs,Sobs)),
             time="time",
             t0=0,
             rprocess=euler.sim(
               step.fun="budmoth_map",
               delta.t=1,
               PACKAGE="pomp.devel"
               ),
             dprocess=onestep.dens(
               dens.fun="budmoth_density",
               PACKAGE="pomp.devel"
               ),
             rmeasure="budmoth_rmeasure",
             dmeasure="budmoth_dmeasure",
             skeleton.type="map",
             skeleton="budmoth_skeleton",
             PACKAGE="pomp.devel",
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
             )
  coef(po) <- unlist(params[d,])
  budmoth.sim[[d]] <- po
}

save(budmoth.sim,file="budmoth.sim.rda",compress="xz")
