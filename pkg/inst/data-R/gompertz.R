require(pomp)

simulate(
         pomp(
              data=data.frame(time=seq(0,100,by=1),Y=NA),
              times="time",
              t0=0,
              rprocess=discrete.time.sim(
                step.fun="_gompertz_simulator"
                ),
              rmeasure="_gompertz_normal_rmeasure",
              dmeasure="_gompertz_normal_dmeasure",
              skeleton.map="_gompertz_skeleton",
              paramnames=c("log.r","log.K","log.sigma","log.tau"),
              statenames=c("X"),
              obsnames=c("Y")
              ),
         params=c(
           log.K=log(1),
           log.r=log(0.1),
           log.sigma=log(0.1),
           log.tau=log(0.1),
           X.0=1
           ),
         nsim=1,
         seed=73691676L
         ) -> gompertz
