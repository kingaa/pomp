library(pomp)

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
           obsnames="reports",
           statenames=c("S","I","R","N","cases"),
           paramnames=c(
             "gamma","mu","iota",
             "beta1","nu",
             "nbasis","degree","period"
             ),
           zeronames=c("cases"),
           measurement.model=reports~binom(size=cases,prob=0.1),
           to.log.transform=c(
             "gamma","nu","mu","iota",
             "beta1","beta2","beta3",
             "rho",
             "S.0","I.0","R.0"
             ),
           parameter.transform=function (params, to.log.transform, ...) {
             params[to.log.transform] <- log(params[to.log.transform])
             params
           },
           parameter.inv.transform=function (params, to.log.transform, ...) {
             params[to.log.transform] <- exp(params[to.log.transform])
             params
           },
           initializer=function(params, t0, ...){
             comp.names <- c("S","I","R")
             icnames <- paste(comp.names,"0",sep=".")
             snames <- c("S","I","R","N","cases")
             fracs <- exp(params[icnames])
             x0 <- numeric(length(snames))
             names(x0) <- snames
             x0["N"] <- params["pop"]
             x0[comp.names] <- round(params["pop"]*fracs/sum(fracs))
             x0
           }
           )

coef(po,transform=TRUE) <- c(
          gamma=24,
          nu=1/70,
          mu=1/70,
          iota=0.1,
          beta1=330,
          beta2=410,
          beta3=490,
          rho=0.1,
          S.0=0.05,
          I.0=1e-4,
          R.0=0.95,
          pop=1000000,
          nbasis=3,
          degree=3,
          period=1
          )


simulate(
         po,
         nsim=1,
         seed=1165270654L
         ) -> gillespie.sir

save(gillespie.sir,file="gillespie.sir.rda",compress=TRUE)
