## This file constructs a pomp object for continuous time Markov chain models
## for pertussis, simulates from it, and plots the data, simulated observations,
## and states.  Simulations from the process model are approximated using
## an Euler approximation.

require(pompExamples)

paramcsv <- "
,birthrate,deathrate,mean.beta,ampl.beta,imports,sigma,gamma,alpha,alpha.ratio,report.prob,boost.prob,polar.prob,foi.mod,popsize,noise.sigma,tau,S.0,E.0,I.0,R1.0,R2.0,seed
SEIR.small,0.02,0.02,450,0.15,10,46,26,0,1,0.3,0,0,0,500000,0,0.01,0.0574148031949802,0.0004081763321755,0.00067028956509212,0.941506730907752,0,1831650124
SEIR.big,0.02,0.02,450,0.15,10,46,26,0,1,0.3,0,0,0,5000000,0,0.01,0.0515635231482973,0.000437143470487014,0.000734641109212043,0.947264692272004,0,908022490
SEIRS.small,0.02,0.02,150,0.15,10,46,26,0.1,1,0.1,0,0,0,500000,0,0.01,0.157360395940609,0.000837874318852172,0.00124181372794081,0.45913512973054,0.381424786282058,1111340406
SEIRS.big,0.02,0.02,150,0.15,10,46,26,0.1,1,0.1,0,0,0,5000000,0,0.01,0.157398354546347,0.00132093662562661,0.0022558671035406,0.457185201591761,0.381839640132725,1751228386
SEIRR.small,0.02,0.02,150,0.15,10,46,26,0.1,1,0.11,0.75,0,0.5,500000,0,0.01,0.128943112158304,0.00068688724266688,0.00114414648269803,0.638074319602244,0.231151534514087,350421545
SEIRR.big,0.02,0.02,150,0.15,10,46,26,0.1,1,0.11,0.75,0,0.5,5000000,0,0.01,0.127128689912424,0.00126497004491763,0.00216092385991776,0.639879739889535,0.229565676293206,748454784
full.small,0.02,0.02,150,0.15,10,46,26,0.1,1,0.1,0.75,0.1,0.5,500000,0.01,0.01,0.132553922599906,0.0010539075727066,0.00166100642162314,0.641737544956371,0.222993618449393,581894515
full.big,0.02,0.02,150,0.15,10,46,26,0.1,1,0.1,0.75,0.1,0.5,5000000,0.01,0.01,0.130980596244438,0.00115096693013597,0.0018994251960431,0.643957103848235,0.222011907781148,301057392
"

tc <- textConnection(paramcsv)
params <- read.csv(tc,row.names=1)
close(tc)

datasets <- rownames(params)
seeds <- params["seed"]
params <- subset(params,select=-seed)

pertussis.sim <- list()
for (d in datasets) {
  po <- pomp(
             data=data.frame(time=seq(from=0,to=20,by=1/52),reports=NA),
             times="time",
             t0=-1/52,
             rprocess = euler.sim(
               step.fun="pertussis_sveirr_EM",
               delta.t=1/52/7,          # Euler stepsize
               PACKAGE="pompExamples"
               ),
             skeleton.type="vectorfield",
             skeleton="pertussis_sveirr_skel",
             PACKAGE="pompExamples",
             paramnames=c(
               "birthrate","deathrate","mean.beta","ampl.beta",
               "imports","sigma","gamma","alpha","alpha.ratio",
               "report.prob","boost.prob","polar.prob",
               "foi.mod","noise.sigma","popsize","tau",
               "S.0","E.0","I.0","R1.0","R2.0"
               ),
             statenames=c("S","E","I","R1","R2","cases","W","err","simpop"),
             zeronames=c("cases","err"),
             ivps=c("S.0","E.0","I.0","R1.0","R2.0"),
             comp.names=c("S","E","I","R1","R2"),
             rmeasure = "negbin_rmeasure",
             dmeasure = "negbin_dmeasure",
             parameter.inv.transform="pertussis_par_untrans",
             parameter.transform="pertussis_par_trans",
             initializer = function (params, t0, statenames, comp.names, ivps, ...) {
               states <- numeric(length(statenames))
               names(states) <- statenames
               ## translate fractions into initial conditions
               frac <- params[ivps]
               states[comp.names] <- round(params['popsize']*frac/sum(frac))
               states["simpop"] <- params["popsize"]
               states
             }
             )
  coef(po) <- unlist(params[d,])
  pertussis.sim[[d]] <- simulate(po,seed=seeds[d,])
}

save(pertussis.sim,file="pertussis.sim.rda",compress="xz")
