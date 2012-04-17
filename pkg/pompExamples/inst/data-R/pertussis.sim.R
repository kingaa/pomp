## This file constructs a pomp object for continuous time Markov chain models
## for pertussis, simulates from it, and plots the data, simulated observations,
## and states.  Simulations from the process model are approximated using
## an Euler approximation.

require(pomp.devel)

params <- read.csv("pertussis-params.csv",sep=";")
datasets <- with(params,paste(model,pop,sep="."))
params <- params[!(names(params)%in%c("model","pop","seed"))]
rownames(params) <- datasets

load("pertussis-simdata.rda")

pertussis.sim <- list()
for (d in datasets) {
  po <- pomp(
             data=subset(pertussis.sim.data,dataset==d,select=c(time,reports)),
             times="time",
             t0=-1/52,
             rprocess = euler.sim(
               step.fun="pertussis_sveirr_EM",
               delta.t=1/52/7,      # Euler stepsize
               PACKAGE="pomp.devel"
               ),
             skeleton.type="vectorfield",
             skeleton="pertussis_sveirr_skel",
             PACKAGE="pomp.devel", # name of the dynamically loadable library where the C functions are
             paramnames=c(
               "birthrate","deathrate","mean.beta","ampl.beta",
               "imports","sigma","gamma","alpha.1","alpha.2","alpha.ratio",
               "report.prob","boost.prob","polar.prob","vacc.prob",
               "foi.mod","noise.sigma","popsize","tau"
               ),
             statenames=c("S","E","I","R1","R2","V","cases","W","err","simpop"),
             zeronames=c("cases","err"),
             ivps=c("S.0","E.0","I.0","R1.0","R2.0","V.0"),
             comp.names=c("S","E","I","R1","R2","V"),
             rmeasure = "negbin_rmeasure",
             dmeasure = "negbin_dmeasure",
             logitvar=c("report.prob","boost.prob","polar.prob","vacc.prob"),
             logvar=c(
               "birthrate","deathrate","mean.beta","ampl.beta","imports","sigma","gamma",
               "alpha.1","alpha.2","alpha.ratio","foi.mod","noise.sigma","tau",
               "S.0","E.0","I.0","R1.0","R2.0","V.0"
               ),
             initializer = function (params, t0, statenames, comp.names, ivps, ...) {
               states <- numeric(length(statenames))
               names(states) <- statenames
               ## translate fractions into initial conditions
               frac <- params[ivps]
               states[comp.names] <- round(params['popsize']*frac/sum(frac))
               states["simpop"] <- params["popsize"]
               states
             },
             parameter.inv.transform=function (params, logitvar, logvar, ivps, ...) {
               params[ivps] <- params[ivps]/sum(params[ivps],na.rm=TRUE)
               params[logitvar] <- qlogis(params[logitvar])
               params[logvar] <- log(params[logvar])
               params
             },
             parameter.transform=function (params, logitvar, logvar, ivps, ...) {
               params[logitvar] <- plogis(params[logitvar])
               params[logvar] <- exp(params[logvar])
               params[ivps] <- params[ivps]/sum(params[ivps],na.rm=TRUE)
               params
             }
             )
  coef(po) <- unlist(params[d,])
  pertussis.sim[[d]] <- po
}

save(pertussis.sim,file="pertussis.sim.rda",compress="xz")
