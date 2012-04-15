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
             "nbasis","degree","period",
             "S.0","I.0","R.0"
             ),
           zeronames=c("cases"),
           comp.names=c("S","I","R"),
           ic.names=c("S.0","I.0","R.0"),
           parameter.transform="_sir_par_trans",
           parameter.inv.transform="_sir_par_untrans",
           initializer=function(params, t0, comp.names, ic.names, ...) {
             snames <- c("S","I","R","cases","W")
             fracs <- params[ic.names]
             x0 <- numeric(length(snames))
             names(x0) <- snames
             x0[comp.names] <- round(params['pop']*fracs/sum(fracs))
             x0
           }
           )

coef(po) <- c(
              gamma=26,mu=0.02,iota=0.01,
              nbasis=3,degree=3,period=1,
              beta1=400,beta2=480,beta3=320,
              beta.sd=1e-3,
              pop=2.1e6,
              rho=0.6,
              S.0=26/400,I.0=0.001,R.0=1-26/400
              )

simulate(po,nsim=1,seed=329343545L) -> euler.sir

save(euler.sir,file="euler.sir.rda",compress="xz")

time(po) <- seq(from=0,to=10,by=1/52)

po <- pomp(
           po,
           statenames=c("S","I","R","N","cases"),
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
           measurement.model=reports~binom(size=cases,prob=rho),
           comp.names=c("S","I","R"),
           ic.names=c("S.0","I.0","R.0"),
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
              nbasis=3,degree=3,period=1,
              beta.sd=0
              )

simulate(
         po,
         nsim=1,
         seed=1165270654L
         ) -> gillespie.sir

save(gillespie.sir,file="gillespie.sir.rda",compress="xz")

tc <- textConnection("
day;reports
1;3
2;8
3;28
4;76
5;222
6;293
7;257
8;237
9;192
10;126
11;70
12;28
13;12
14;5
")

flu <- read.csv2(file=tc)
close(tc)

po <- pomp(
           data=flu,
           times="day",
           t0=0,
           rprocess=euler.sim(
             step.fun="_sir_euler_simulator",
             delta.t=1/12,
             PACKAGE="pomp"
             ),
           skeleton.type="vectorfield",
           skeleton="_sir_ODE",
           measurement.model=reports~norm(mean=rho*cases,sd=1+sigma*cases),
           PACKAGE="pomp",
           obsnames = c("reports"),
           statenames=c("S","I","R","cases","W"),
           paramnames=c(
             "gamma","mu","iota",
             "beta","beta.sd","pop","rho",
             "nbasis","degree","period",
             "S.0","I.0","R.0"
             ),
           zeronames=c("cases"),
           comp.names=c("S","I","R"),
           ic.names=c("S.0","I.0","R.0"),
           logvar=c(
             "beta","gamma","mu","iota","sigma","beta.sd",
             "S.0","I.0","R.0"
             ),
           logitvar="rho",
           parameter.inv.transform=function (params, logvar, logitvar, ...) {
             params[logvar] <- log(params[logvar])
             params[logitvar] <- qlogis(params[logitvar])
             params
           },
           parameter.transform=function (params, logvar, logitvar, ...) {
             params[logvar] <- exp(params[logvar])
             params[logitvar] <- plogis(params[logitvar])
             params
           },
           initializer=function(params, t0, comp.names, ic.names, ...) {
             snames <- c("S","I","R","cases","W")
             fracs <- params[ic.names]
             x0 <- numeric(length(snames))
             names(x0) <- snames
             x0[comp.names] <- round(params['pop']*fracs/sum(fracs))
             x0
           }
           )

coef(po) <- c(
              gamma=1/3,mu=0.0,iota=0.0,
              nbasis=1,degree=0,period=1,
              beta=1.4,
              beta.sd=0,
              pop=1400,
              rho=0.9,sigma=3.6,
              S.0=0.999,I.0=0.001,R.0=0
              )

bbs <- po
save(bbs,file="bbs.rda",compress="xz")
