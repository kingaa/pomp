library(pomp)

params <- c(
            nu=1/70,
            mu=1/70,
            beta1=330,
            beta2=410,
            beta3=490,
            gamma=24,
            iota=0.1,
            S.0=0.05,
            I.0=1e-4,
            R.0=0.95,
            N.0=1000000
            )

seasonality <- data.frame(
                          time=seq(0,10,by=1/52/10),
                          seas=periodic.bspline.basis(seq(0,10,by=1/52/10),nbasis=3,degree=3,period=1)
                          )
simulate(
         pomp(
              data=data.frame(
                time=seq(from=0,to=2,by=1/52),
                reports=NA
                ),
              times="time",
              t0=0,
              covar=seasonality,
              tcovar="time",
              rprocess=gillespie.sim(
                rate.fun=function(j, x, t, params, covars, ...) {
                  switch(
                         j,
                         params["nu"]*x["N"], # birth
                         params["mu"]*x["S"], # susceptible death
                         {                         # infection
                           beta <- params[c("beta1","beta2","beta3")]%*%covars[c("seas.1","seas.2","seas.3")]
                           (beta*x["I"]+params["iota"])*x["S"]/x["N"]
                         },
                         params["mu"]*x["I"], # infected death
                         params["gamma"]*x["I"], # recovery
                         params["mu"]*x["R"], # recovered death
                         stop("unrecognized event ",j)
                         )
                },
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
              zeronames=c("cases"),
              measurement.model=reports~binom(size=cases,prob=0.1),
              initializer=function(params, t0, ...){
                comp.names <- c("S","I","R")
                icnames <- paste(comp.names,"0",sep=".")
                snames <- c("S","I","R","N","cases")
                fracs <- params[icnames]
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
         ) -> gsir


tail(as(gsir,"data.frame"))

data(gillespie.sir)

tail(as(simulate(gillespie.sir,times=time(gsir),t0=timezero(gsir),seed=1165270654L),"data.frame"))
