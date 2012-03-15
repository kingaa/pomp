require(pomp)

po <- pomp(
           data=data.frame(
             time=seq(from=1/52,to=4,by=1/52),
             reports=NA
             ),
           times="time",
           t0=0,
           rprocess=euler.sim(
             step.fun="_sir_euler_simulator",
             delta.t=1/52/20,
             PACKAGE="pomp"
             ),
           skeleton.type="vectorfield",
           skeleton="_sir_ODE",
           rmeasure="_sir_binom_rmeasure",
           dmeasure="_sir_binom_dmeasure",
           PACKAGE="pomp",
           obsnames = c("reports"),
           statenames=c("S","I","R","cases","W"),
           paramnames=c(
             "gamma","mu","iota",
             "beta1","beta.sd","pop","rho",
             "nbasis","degree","period"
             ),
           zeronames=c("cases"),
           comp.names=c("S","I","R"),
           to.log.transform=c(
             "gamma","mu","iota",
             "beta1","beta2","beta3","beta.sd",
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
           initializer=function(params, t0, comp.names, ...) {
             ic.names <- paste(comp.names,".0",sep="")
             snames <- c("S","I","R","cases","W")
             fracs <- exp(params[ic.names])
             x0 <- numeric(length(snames))
             names(x0) <- snames
             x0[comp.names] <- round(params['pop']*fracs/sum(fracs))
             x0["cases"] <- 0
             x0
           }
           )
coef(po,transform=TRUE) <- c(gamma=26,mu=0.02,iota=0.01,
          nbasis=3,degree=3,period=1,
          beta1=1200,beta2=1800,beta3=600,
          beta.sd=1e-3,
          pop=2.1e6,
          rho=0.6,
          S.0=26/1200,I.0=0.001,R.0=1-0.001-26/1200
          )

simulate(po,nsim=1,seed=329348545L) -> euler.sir

save(euler.sir,file="euler.sir.rda",compress="xz")

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
           measurement.model=reports~binom(size=cases,prob=exp(rho)),
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

save(gillespie.sir,file="gillespie.sir.rda",compress="xz")
