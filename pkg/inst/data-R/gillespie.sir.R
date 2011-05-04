library(pomp)

params <- c(
            nu=log(1/70),
            mu=log(1/70),
            beta1=log(330),
            beta2=log(410),
            beta3=log(490),
            gamma=log(24),
            iota=log(0.1),
            rho=log(0.1),
            S.0=log(0.05),
            I.0=log(1e-4),
            R.0=log(0.95),
            N.0=1000000,
            cases.0=0,
            nbasis=3,
            degree=3,
            period=1
            )

simulate(
         pomp(
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
              initializer=function(params, t0, ...){
                comp.names <- c("S","I","R")
                icnames <- paste(comp.names,"0",sep=".")
                snames <- c("S","I","R","N","cases")
                fracs <- exp(params[icnames])
                x0 <- numeric(length(snames))
                names(x0) <- snames
                x0["N"] <- params["N.0"]
                x0[comp.names] <- round(params['N.0']*fracs/sum(fracs))
                x0
              }
              ),
         params=params,
         nsim=1,
         seed=1165270654L
         ) -> gillespie.sir

save(gillespie.sir,file="gillespie.sir.rda")
