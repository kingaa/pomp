## blowfly model, with general dt
## here, set up for dt=1 and dt=2
## dt is hard-coded, and initial values are customized for each dt

require(pomp)

## following xia and tong, the delay is treated as fixed at 14 days
## xia and tong claim to be using tau=8 bidays, but on inspection 
## their Euler method is really tau=7 bidays

raw.data <- subset(
                   read.csv2("blowflies.csv",comment.char="#"),
                   set==4
                   )

pomp(
     data=subset(raw.data[c("day","y")],day>14&day<400),
     times="day",
     t0=14,
     rprocess=discrete.time.sim(
       step.fun="_blowfly_model_simulator",
       delta.t=1,
       ),
     paramnames=c("log.P","log.N0","log.delta","log.sigma.P","log.sigma.d","tau","log.sigma.y"),
     statenames=c("N1","R","S","e","eps"),
     obsnames=c("y"),
     measurement.model=y~nbinom(mu=N1,size=exp(-2*log.sigma.y)),
     y.init=with( ## initial data
       raw.data,
       approx(
              x=day,
              y=y,
              xout=seq(from=0,to=14,by=1),
              rule=2
              )$y
       ),
#     y.init=c(948, 948, 942, 930, 911, 885, 858, 833.7, 801, 748.3, 676, 589.8, 504, 434.9, 397),
     initializer=function (params, t0, y.init, ...) {
       ntau <- length(y.init)
       n <- y.init[ntau:1]
       names(n) <- paste("N",seq_len(ntau),sep="")
       c(n,R=0,S=0,e=0,eps=0)
     }
     ) -> blowflies1

pomp(
     data=subset(raw.data[c("day","y")],day>14&day<400),
     times="day",
     t0=14,
     rprocess=discrete.time.sim(
       step.fun="_blowfly_model_simulator",
       delta.t=2,
       ),
     paramnames=c("log.P","log.N0","log.delta","log.sigma.P","log.sigma.d","tau","log.sigma.y"),
     statenames=c("N1","R","S","e","eps"),
     obsnames=c("y"),
     measurement.model=y~nbinom(mu=N1,size=exp(-2*log.sigma.y)),
     y.init=with( ## initial data
       raw.data,
       approx(
              x=day,
              y=y,
              xout=seq(from=0,to=14,by=2),
              rule=2
              )$y
       ),
     #y.init=c(948, 942, 911, 858, 801, 676, 504, 397),
     initializer=function (params, t0, y.init, ...) {
       ntau <- length(y.init)
       n <- y.init[ntau:1]
       names(n) <- paste("N",seq_len(ntau),sep="")
       c(n,R=0,S=0,e=0,eps=0)
     }
     ) -> blowflies2

## mle from search to date
coef(blowflies1) <- c(
                    log.P = 1.189 , 
                    log.delta = -1.828 , 
                    log.N0 = 6.522 , 
                    log.sigma.P = 0.301 , 
                    log.sigma.d = -0.292 , 
                    log.sigma.y = -3.625 , 
                    tau = 14 
                    )

## mle from search to date
coef(blowflies2) <- c(
                    log.P = 1.005 , 
                    log.delta = -1.75 , 
                    log.N0 = 6.685 , 
                    log.sigma.P = 0.366 , 
                    log.sigma.d = -0.274 , 
                    log.sigma.y = -4.524 , 
                    tau = 7 
                    )

test <- FALSE
if(test){
  sim1 <- simulate(blowflies1,nsim=1)
  plot(obs(sim1)['y',],ty='l')
  lines(obs(blowflies1)['y',],lty="dashed")
  states(sim1)[,1]

  sim2 <- simulate(blowflies2,nsim=1)
  plot(obs(sim2)['y',],ty='l')
  lines(obs(blowflies2)['y',],lty="dashed")
  states(sim2)[,1]

  ## check that it matches the deterministic skeleton when noise is small
  params.1.skel <- coef(blowflies1)
  params.1.skel["log.sigma.P"] <- log(0.00001)
  params.1.skel["log.sigma.d"] <- log(0.00001)
  params.1.skel["log.sigma.y"] <- log(0.00001)
  simulate(blowflies1,params=params.1.skel,nsim=1,seed=73691676L) -> b1.skel
  plot(obs(blowflies1)['y',],ty='l',lty="dashed")
  lines(obs(b1.skel)['y',],ty='l')
  
} 

save(blowflies1,blowflies2,file="blowflies.rda",compress="xz")
