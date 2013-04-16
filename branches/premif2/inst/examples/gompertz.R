require(pomp)

po <- pomp(
           data=data.frame(time=seq(0,100,by=1),Y=NA),
           times="time",
           t0=0,
           rprocess=discrete.time.sim(
             step.fun="_gompertz_simulator"
             ),
           rmeasure="_gompertz_normal_rmeasure",
           dmeasure="_gompertz_normal_dmeasure",
           skeleton.type="map",
           skeleton="_gompertz_skeleton",
           paramnames=c("r","K","sigma","tau"),
           statenames=c("X"),
           obsnames=c("Y"),
           parameter.transform=function(params,...){
             exp(params)
           },
           parameter.inv.transform=function(params,...){
             log(params)
           }
           )

coef(po) <- c(K=1,r=0.1,sigma=0.1,tau=0.1,X.0=1)

simulate(po,nsim=1,seed=299438676L) -> gompertz

assign("gompertz",gompertz,envir=.GlobalEnv)
c("gompertz")
