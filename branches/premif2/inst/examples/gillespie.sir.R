require(pomp)

po <- pomp(
           data=data.frame(
             time=seq(from=0,to=10,by=1/52),
             reports=NA
             ),
           times="time",
           t0=0,
           rprocess=gillespie.sim(
             rate.fun="_sir_rates",
             PACKAGE="pomp",
             v=cbind(
               birth=c(1,0,0,1,0),
               sdeath=c(-1,0,0,-1,0),
               infection=c(-1,1,0,0,0),
               ideath=c(0,-1,0,-1,0),
               recovery=c(0,-1,1,0,1),
               rdeath=c(0,0,-1,-1,0)
               ),
             d=cbind(
               birth=c(0,0,0,1,0),
               sdeath=c(1,0,0,0,0),
               infection=c(1,1,0,1,0),
               ideath=c(0,1,0,0,0),
               recovery=c(0,1,0,0,0),
               rdeath=c(0,0,1,0,0)
               )
             ),
           skeleton.type="vectorfield",
           skeleton="_sir_ODE",
           measurement.model=reports~binom(size=cases,prob=rho),
           PACKAGE="pomp",
           obsnames = c("reports"),
           statenames=c("S","I","R","N","cases"),
           paramnames=c(
             "gamma","mu","iota",
             "beta1","beta.sd","pop","rho",
             "S.0","I.0","R.0"
             ),
           zeronames=c("cases"),
           comp.names=c("S","I","R"),
           ic.names=c("S.0","I.0","R.0"),
           parameter.transform="_sir_par_trans",
           parameter.inv.transform="_sir_par_untrans",
           nbasis=3L,
           degree=3L,
           period=1.0,
           initializer=function(params, t0, comp.names, ic.names, ...) {
             x0 <- numeric(5)
             names(x0) <- c("S","I","R","N","cases")
             fracs <- params[ic.names]
             x0["N"] <- params["pop"]
             x0[comp.names] <- round(params["pop"]*fracs/sum(fracs))
             x0
           }
           )

coef(po) <- c(
              gamma=24,mu=1/70,iota=0.1,
              beta1=330,beta2=410,beta3=490,
              rho=0.1,
              S.0=0.05,I.0=1e-4,R.0=0.95,
              pop=1000000,
              beta.sd=0
              )

simulate(
         po,
         nsim=1,
         seed=1165270654L
         ) -> gillespie.sir

assign("gillespie.sir",gillespie.sir,envir=.GlobalEnv)
c("gillespie.sir")
