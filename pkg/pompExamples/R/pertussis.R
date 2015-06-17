## This file constructs a pomp object for one of several continuous time Markov chain models
## for pertussis containing simulated data.
## Simulations from the process model are approximated using an Euler approximation.

pertussis.sim <- function (which) {
  if (missing(which)) {
    datasets <- c(
                  "SEIR.small","SEIR.big",
                  "SEIRS.small","SEIRS.big",
                  "SEIRR.small","SEIRR.big",
                  "full.small","full.big"
                  )

    cat("available datasets:",sQuote(datasets),"\n")
    invisible(datasets)
  } else {
    which <- as.character(substitute(which))
    suppressMessages(
                     simulate(
                              pomp(
                                   data=data.frame(time=seq(from=0,to=20,by=1/52),reports=NA),
                                   times="time",
                                   t0=-1/52,
                                   params=switch(
                                     which,
                                     SEIR.small=c(
                                       birthrate=0.02, deathrate=0.02, mean.beta=450, ampl.beta=0.15,
                                       imports=10, sigma=46, gamma=26, alpha=0, alpha.ratio=1,
                                       report.prob=0.3, boost.prob=0, polar.prob=0, foi.mod=0,
                                       popsize=5e+5, noise.sigma=0, tau=0.01,
                                       S.0=0.0574148031949802, E.0=0.0004081763321755, I.0=0.00067028956509212,
                                       R1.0=0.941506730907752, R2.0=0),
                                     SEIR.big=c(
                                       birthrate=0.02, deathrate=0.02, mean.beta=450, ampl.beta=0.15,
                                       imports=10, sigma=46, gamma=26, alpha=0, alpha.ratio=1,
                                       report.prob=0.3, boost.prob=0, polar.prob=0, foi.mod=0,
                                       popsize=5e+6, noise.sigma=0, tau=0.01,
                                       S.0=0.0515635231482973, E.0=0.000437143470487014, I.0=0.000734641109212043,
                                       R1.0=0.947264692272004, R2.0=0),
                                     SEIRS.small=c(
                                       birthrate=0.02, deathrate=0.02, mean.beta=150, ampl.beta=0.15,
                                       imports=10, sigma=46, gamma=26, alpha=0.1, alpha.ratio=1,
                                       report.prob=0.1, boost.prob=0, polar.prob=0, foi.mod=0,
                                       popsize=5e+5, noise.sigma=0, tau=0.01,
                                       S.0=0.157360395940609, E.0=0.000837874318852172, I.0=0.00124181372794081,
                                       R1.0=0.45913512973054, R2.0=0.381424786282058),
                                     SEIRS.big=c(
                                       birthrate=0.02, deathrate=0.02, mean.beta=150, ampl.beta=0.15,
                                       imports=10, sigma=46, gamma=26, alpha=0.1, alpha.ratio=1,
                                       report.prob=0.1, boost.prob=0, polar.prob=0, foi.mod=0,
                                       popsize=5e+6, noise.sigma=0, tau=0.01,
                                       S.0=0.157398354546347, E.0=0.00132093662562661, I.0=0.0022558671035406,
                                       R1.0=0.457185201591761, R2.0=0.381839640132725),
                                     SEIRR.small=c(
                                       birthrate=0.02, deathrate=0.02, mean.beta=150, ampl.beta=0.15,
                                       imports=10, sigma=46, gamma=26, alpha=0.1, alpha.ratio=1,
                                       report.prob=0.11, boost.prob=0.75, polar.prob=0, foi.mod=0.5,
                                       popsize=5e+5, noise.sigma=0, tau=0.01,
                                       S.0=0.128943112158304, E.0=0.00068688724266688, I.0=0.00114414648269803,
                                       R1.0=0.638074319602244, R2.0=0.231151534514087),
                                     SEIRR.big=c(
                                       birthrate=0.02, deathrate=0.02, mean.beta=150, ampl.beta=0.15,
                                       imports=10, sigma=46, gamma=26, alpha=0.1, alpha.ratio=1,
                                       report.prob=0.11, boost.prob=0.75, polar.prob=0, foi.mod=0.5,
                                       popsize=5e+6, noise.sigma=0, tau=0.01,
                                       S.0=0.127128689912424, E.0=0.00126497004491763, I.0=0.00216092385991776,
                                       R1.0=0.639879739889535, R2.0=0.229565676293206),
                                     full.small=c(
                                       birthrate=0.02, deathrate=0.02, mean.beta=150, ampl.beta=0.15,
                                       imports=10, sigma=46, gamma=26, alpha=0.1, alpha.ratio=1,
                                       report.prob=0.1, boost.prob=0.75, polar.prob=0.1, foi.mod=0.5,
                                       popsize=5e+5, noise.sigma=0.01, tau=0.01,
                                       S.0=0.132553922599906, E.0=0.0010539075727066, I.0=0.00166100642162314,
                                       R1.0=0.641737544956371, R2.0=0.222993618449393),
                                     full.big=c(
                                       birthrate=0.02, deathrate=0.02, mean.beta=150, ampl.beta=0.15,
                                       imports=10, sigma=46, gamma=26, alpha=0.1, alpha.ratio=1,
                                       report.prob=0.1, boost.prob=0.75, polar.prob=0.1, foi.mod=0.5,
                                       popsize=5e+6, noise.sigma=0.01, tau=0.01,
                                       S.0=0.130980596244438, E.0=0.00115096693013597, I.0=0.0018994251960431,
                                       R1.0=0.643957103848235, R2.0=0.222011907781148),
                                     stop("unrecognized dataset ",sQuote(which),call.=FALSE)
                                     ),
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
                                   toEstimationScale="pertussis_par_untrans",
                                   fromEstimationScale="pertussis_par_trans",
                                   varnames=c("S","E","I","R1","R2","cases","W","err","simpop"),
                                   initializer = function (params, t0, varnames, comp.names, ivps, ...) {
                                     states <- numeric(length(varnames))
                                     names(states) <- varnames
                                     ## translate fractions into initial conditions
                                     frac <- params[ivps]
                                     states[comp.names] <- round(params['popsize']*frac/sum(frac))
                                     states["simpop"] <- params["popsize"]
                                     states
                                   }
                                   ),
                              seed=switch(
                                which,
                                SEIR.small=1831650124L,
                                SEIR.big=908022490L,
                                SEIRS.small=1111340406L,
                                SEIRS.big=1751228386L,
                                SEIRR.small=350421545L,
                                SEIRR.big=748454784L,
                                full.small=581894515L,
                                full.big=301057392L,
                                stop("unrecognized dataset ",sQuote(which),call.=FALSE)
                                )
                              )
                     )
  }
}
