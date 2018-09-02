options(digits=3)
png(filename="spect-%02d.png",res=100)

library(pomp)

pompExample(gompertz)
set.seed(362083261L)

sp <- spect(gompertz,kernel.width=3,nsim=100,seed=838775L)
summary(sp)
plot(sp)

spect(sp,kernel.width=5,seed=838775L) -> sp2
plot(sp2)

spect(sp,detrend="mean") -> sp3
spect(sp,detrend="linear") -> sp4
spect(sp,detrend="quadratic") -> sp5

theta <- as.list(coef(sp))
theta$r <- 25
spect(sp,params=theta) -> sp4
plot(sp4,quantiles=c(0.5,0.95))
plot(sp4,plot.data=FALSE,quantile.styles=list(col=1:5))
plot(sp4,plot.data=FALSE,quantile.styles=list(col="grey10",lty=1:5))
plot(sp4,plot.data=FALSE,quantile.styles=list(col="grey10",lty=1:3))
try(plot(sp4,plot.data=FALSE,quantile.styles=c(col="grey10",lty=1:3)))
plot(sp4,plot.data=TRUE,data.styles=list(col="red",lty=1))
try(plot(sp4,plot.data=TRUE,data.styles=c(col="red",lty=1)))

try(spect())
try(spect("bob"))

try(spect(sp,kernel.width=-3))
try(spect(sp,kernel.width=NA))
try(spect(sp,kernel.width=NULL))

try(spect(sp,nsim=-100))
try(spect(sp,nsim=Inf))
try(spect(sp,nsim=NULL))

sp4@data[17] <- NA
try(spect(sp4))

try(spect(sp3,rmeasure=function(t,X,...){
  if (t==13) c(Y=NA) else c(Y=X)
}))

time(sp3) <- c(0:7,10:40)
try(spect(sp3))

try(spect(sp2,rmeasure=function(...) stop("yikes!")))

simulate(times=1:100,t0=0,
  rprocess=euler(Csnippet("
       x = rnorm(0,1);
       y = rnorm(0,1);
       z = rnorm(0,1);"),
    delta.t=1),
  rmeasure=Csnippet("
       a = rnorm(x,1);
       b = rnorm(y,1);
       c = rnorm(z,1);"),
  rinit=Csnippet("x = y = z = 0;"),
  obsnames=c("a","b","c"),
  statenames=c("x","y","z"),
  params=c()) -> bob

plot(spect(bob,kernel.width=3,nsim=500),
  data.style=list(lwd=c(2,3),lty=2,col='red'))

dev.off()
