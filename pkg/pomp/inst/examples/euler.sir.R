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
              beta1=400,beta2=480,beta3=320,
              beta.sd=1e-3,
              pop=2.1e6,
              rho=0.6,
              S.0=26/400,I.0=0.001,R.0=1-26/400
              )

simulate(po,nsim=1,seed=329343545L) -> euler.sir

assign("euler.sir",euler.sir,envir=.GlobalEnv)
c("euler.sir")
