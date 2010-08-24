require(pomp)

simulate(
         pomp(
              data=data.frame(time=seq(0,50,by=1),y=NA),
              times="time",
              t0=0,
              rprocess=discrete.time.sim(
                step.fun="ricker_simulator"
                ),
              rmeasure="poisson_rmeasure",
              dmeasure="poisson_dmeasure",
              skeleton.map="ricker_skeleton",
              paramnames=c("log.r","log.sigma","log.phi"),
              statenames=c("N","e"),
              obsnames=c("y")
              ),
         params=c(
           log.r=3.8,
           log.sigma=log(0.3),
           log.phi=log(10),
           N.0=7,
           e.0=0
           ),
         nsim=1,
         seed=73691676L
         ) -> ricker
