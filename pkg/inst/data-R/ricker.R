require(pomp)

simulate(
         pomp(
              data=data.frame(time=seq(0,50,by=1),y=NA),
              times="time",
              t0=0,
              rprocess=discrete.time.sim(
                step.fun="_ricker_simulator"
                ),
              rmeasure="_ricker_poisson_rmeasure",
              dmeasure="_ricker_poisson_dmeasure",
              skeleton.type="map",
              skeleton="_ricker_skeleton",
              paramnames=c("log.r","sigma","phi"),
              statenames=c("N","e"),
              obsnames=c("y"),
              parameter.inv.transform=function(params,...) {
                params <- c(params["log.r"],log(params[c("sigma","phi","N.0")]),params["e.0"])
                names(params) <- c("log.r","log.sigma","log.phi","log.N.0","e.0")
                params
              },
              parameter.transform=function(params,...) {
                params <- c(params["log.r"],exp(params[c("log.sigma","log.phi","log.N.0")]),params["e.0"])
                names(params) <- c("log.r","sigma","phi","N.0","e.0")
                params
              },
              PACKAGE="pomp"
              ),
         params=c(
           log.r=3.8,
           sigma=0.3,
           phi=10,
           N.0=7,
           e.0=0
           ),
         nsim=1,
         seed=73691676L
         ) -> ricker

save(ricker,file="ricker.rda",compress="xz")
