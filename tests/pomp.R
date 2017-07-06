library(pomp)

png(filename="pomp-%02d.png",res=100)

set.seed(3488585L)

pompExample(ricker)
y1 <- obs(simulate(ricker,seed=1066L))
r2 <- pomp(ricker,measurement.model=y~pois(lambda=N*phi))
coef(r2) <- coef(ricker)
y2 <- obs(simulate(r2,seed=1066L))
max(abs(y1-y2))
r3 <- pomp(
           ricker,
           dmeasure="_ricker_poisson_dmeasure",
           PACKAGE="pomp",
           paramnames=c("r","sigma","phi","c"),
           statenames=c("N","e"),
           obsnames=c("y")
           )
coef(r3) <- coef(r2)
y3 <- obs(simulate(r3,seed=1066L))
max(abs(y1-y3))
r4 <- pomp(
           r2,
           rmeasure="_ricker_poisson_rmeasure",
           PACKAGE="pomp",
           paramnames=c("r","sigma","phi","c"),
           statenames=c("N","e"),
           obsnames=c("y")
           )
coef(r4) <- coef(r2)
y4 <- obs(simulate(r4,seed=1066L))
max(abs(y1-y4))

dat <- as.data.frame(ricker)
try(pomp(dat) -> po)
try(pomp(dat,times="time",t0=0,covar=dat) -> po)
try(pomp(dat,times="time",t0=0,covar=dat,tcovar=3) -> po)
pomp(dat,times="time",t0=0,covar=dat,tcovar=1) -> po
pomp(dat,times=1,t0=0,covar=dat,tcovar=1) -> po
try(pomp(dat,times="time",t0=0,covar=dat,tcovar="bob") -> po)
try(pomp(dat,times="time",t0=0,covar=dat,tcovar=1,covarnames="henry") -> po)
try(pomp(dat,times="time",t0=0,fromEstimationScale=identity) -> po)
try(pomp(dat$y,times=dat$time,t0=0,skeleton.type="map",skelmap.delta.t=1) -> po)
try(pomp(dat,times=1,t0=0,skeleton.type="map",skelmap.delta.t=1) -> po)
plot(po)
pomp(dat,times="time",t0=0,skeleton=function(x,t,params,...){x}) -> po
try(pomp())
try(pomp(as.matrix(dat),times=dat$time,t0=0) -> po)
try(pomp(t(as.matrix(dat)),times=dat$time,t0=0,
         skeleton.type="map",skelmap.delta.t=1) -> po)
try(pomp(dat,times=1,t0=0,skeleton.type="map",skelmap.delta.t=1) -> po)
pomp(dat,times="time",t0=0,skeleton=function(x,t,params,...){x}) -> po
try(pomp(dat,times="bob",t0=0,skeleton=function(x,t,params,...){x}) -> po)
try(pomp(dat,times=11,t0=0,skeleton=function(x,t,params,...){x}) -> po)
try(pomp(dat,times=1.1,t0=0,skeleton=function(x,t,params,...){x}) -> po)
try(pomp(dat$y,times=dat$time,t0=0) -> po)
try(pomp(dat,times=1,t0=0,skeleton.type="map",skelmap.delta.t=1) -> po)
try(pomp(dat,times=1,t0=0,skelmap.delta.t=1) -> po)
pomp(dat,times=1,t0=0,skeleton=function(x,t,params,...){x}) -> po
try(pomp(dat$y,times=dat$time[1:10],t0=0) -> po)
try(pomp(ricker,skeleton.type="map",skelmap.delta.t=1) -> po)
try(pomp(ricker,skelmap.delta.t=1) -> po)
try(pomp(ricker,skeleton=identity(identity)) -> po)
try(pomp(ricker,toEstimationScale=identity) -> po)
try(pomp(ricker,fromEstimationScale=identity) -> po)
try(pomp("banana"))
pomp(ricker,measurement.model=y~pois(N),rmeasure=Csnippet("y=rpois(N);")) -> po

xdat <- ricker@data
nm <- dimnames(xdat)
dim(xdat) <- c(nrow(xdat),1,ncol(xdat))
dimnames(xdat) <- list(nm[[1]],NULL,nm[[2]])
try(pomp(xdat,times=dat$time,t0=0) -> po)

try(pompExample(bob))
try(pompExample("bob"))
pompExample("ricker")
pomp(ricker) -> ricker
pomp(ricker,rmeasure=Csnippet("y=rpois(N);"),statenames="N") -> po
simulate(po) -> po
## force recompile
file.remove(list.files(path=file.path(tempdir(),Sys.getpid()),
                       pattern=paste0(po@solibs[[1]]$name,".*"),
                       full.names=TRUE))
simulate(po) -> po

plot(po,yax.flip=TRUE)

dev.off()
