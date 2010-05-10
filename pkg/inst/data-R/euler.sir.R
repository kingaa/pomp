require(pomp)

simulate(
         pomp(
              times=seq(1/52,4,by=1/52),
              data=rbind(measles=numeric(52*4)),
              t0=0,
              statenames=c("S","I","R","cases","W"),
              paramnames=c(
                "gamma","mu","iota",
                "beta1","beta.sd","pop","rho",
                "nbasis","degree","period"
                ),
              zeronames=c("cases"),
              comp.names=c("S","I","R"),
              rprocess=euler.sim(
                step.fun="sir_euler_simulator",
                delta.t=1/52/20,
                PACKAGE="pomp"
                ),
              skeleton.vectorfield="sir_ODE",
              rmeasure="binom_rmeasure",
              dmeasure="binom_dmeasure",
              PACKAGE="pomp",
              initializer=function(params, t0, comp.names, ...){
                p <- exp(params)
                snames <- c("S","I","R","cases","W")
                fracs <- p[paste(comp.names,"0",sep=".")]
                x0 <- numeric(length(snames))
                names(x0) <- snames
                x0[comp.names] <- round(p['pop']*fracs/sum(fracs))
                x0
              }
              ),
         params=c(
           gamma=log(26),mu=log(0.02),iota=log(0.01),
           nbasis=3,degree=3,period=1, # NB: all parameters are log-transformed but these
           beta1=log(1200),beta2=log(1800),beta3=log(600),
           beta.sd=log(1e-3),
           pop=log(2.1e6),
           rho=log(0.6),
           S.0=log(26/1200),I.0=log(0.001),R.0=log(1-0.001-26/1200)
           ),
         nsim=1,
         seed=329348545L
         ) -> euler.sir
