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
              paramnames=c("r","sigma","phi"),
              statenames=c("N","e"),
              obsnames=c("y"),
              parameter.inv.transform=function(params,...) {
                params[c("r","sigma","phi","N.0")] <- log(params[c("r","sigma","phi","N.0")])
                params
              },
              parameter.transform=function(params,...) {
                params[c("r","sigma","phi","N.0")] <- exp(params[c("r","sigma","phi","N.0")])
                params
              },
              PACKAGE="pomp"
              ),
         params=c(
           r=exp(3.8),
           sigma=0.3,
           phi=10,
           N.0=7,
           e.0=0
           ),
         nsim=1,
         seed=73691676L
         ) -> ricker

assign("ricker",ricker,envir=.GlobalEnv)
c("ricker")
