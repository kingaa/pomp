require(pompExamples)

load("budmoth-params.rda")

datasets <- rownames(params)
seeds <- params["seed"]
params <- subset(params,select=-seed)

budmoth.sim <- list()
for (d in datasets) {
  po <- pomp(
             data=data.frame(
               time=seq(from=0,to=60,by=1),
               Qobs=NA,Nobs=NA,Sobs=NA
               ),
             time="time",
             t0=0,
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
             )
  coef(po) <- unlist(params[d,])
  budmoth.sim[[d]] <- simulate(po,seed=seeds[d,])
}

save(budmoth.sim,file="budmoth.sim.rda",compress="xz")
