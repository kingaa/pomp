require(pomp)

simulate(
         pomp( 
              times=seq(1,101),
              data=rbind(
                y1=rep(0,101),
                y2=rep(0,101)
                ),
              t0=1,
              rprocess=discrete.time.sim("ou2_step",PACKAGE="pomp"),
              dprocess=onestep.dens("ou2_pdf",PACKAGE="pomp"),
              dmeasure = "ou2_dmeasure",
              rmeasure = "ou2_rmeasure",
              skeleton.type="map",
              skeleton = "ou2_skel",
              PACKAGE="pomp",
              paramnames = c(
                "alpha.1","alpha.2","alpha.3","alpha.4",
                "sigma.1","sigma.2","sigma.3",
                "tau"
                ),
              statenames = c("x1","x2"),
              obsnames = c("y1","y2")
              ),
         params=c(
           alpha.1=0.8, alpha.2=-0.5, alpha.3=0.3, alpha.4=0.9,
           sigma.1=3, sigma.2=-0.5, sigma.3=2,
           tau=1, 
           x1.0=-3, x2.0=4
           ),
         nsim=1,
         seed=377456545L
         ) -> ou2

ou2 <- window(ou2,end=100)
timezero(ou2) <- 0

save(ou2,file="ou2.rda")
