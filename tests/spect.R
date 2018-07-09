library(pomp)
set.seed(362083261L)

png(filename="spect-%02d.png",res=100)

pompExample(ou2)

gm1 <- spect.match(ou2,kernel.width=3,detrend="mean",nsim=50,
                   est=c("alpha.1","alpha.4"),reltol=1e-3,maxit=100,
                   method="Nelder-Mead")
gm2 <- spect.match(ou2,kernel.width=3,detrend="mean",nsim=49,
                   est=c("alpha.1","alpha.4"),reltol=1e-3,maxit=100,
                   method="Nelder-Mead")
gm3 <- spect.match(ou2,kernel.width=3,detrend="linear",nsim=50,
                   est=c("alpha.1","alpha.4"),reltol=1e-3,maxit=100,
                   method="Nelder-Mead")
gm4 <- spect.match(ou2,kernel.width=5,detrend="quadratic",nsim=50,
                   est=c("alpha.1","alpha.4"),reltol=1e-3,maxit=100,
                   method="Nelder-Mead")
gm5 <- spect.match(ou2,kernel.width=3,nsim=50,
                   est=c("alpha.1","alpha.4"),reltol=1e-3,maxit=100,
                   method="Nelder-Mead",start=as.list(coef(ou2)))
plot(gm4,y=NA)
plot(gm4,plot.data=FALSE)
plot(gm4,plot.data=FALSE,quantile.styles=list(col=1:5))
plot(gm4,quantile.styles=list(lwd=c(1,2,3)))
try(plot(gm4,plot.data=FALSE,quantile.styles=c(lwd=c(1,2,3))))
try(plot(gm4,data.styles=c(lty=c(1,2,3))))
plot(gm4,data.styles=list(lty=c(1,2,3)))
plot(gm4,data.styles=list(lty=1))

pompExample(ricker)

sp <- spect(ricker,kernel.width=3,nsim=100,seed=838775L)
invisible(summary(sp))
spp <- spect.match(sp,est="")
invisible(summary(spp))

spp <- spect.match(sp,nsim=100,est=c("sigma","phi"),reltol=1e-3,maxit=100)
invisible(summary(spp))

po <- ricker
coef(po,"r") <- 5
sp <- spect(po,kernel.width=1,nsim=100,seed=838775L)
invisible(summary(sp))

po <- ricker
theta <- as.list(coef(po))
theta["phi"] <- 30
sp <- spect(po,kernel.width=11,nsim=100,seed=838775L,params=theta)
invisible(summary(sp))
plot(simulate(sp),variables="y")
sp <- spect(sp)

po <- ricker
try(spect(po))
try(spect(po,kernel.width=3,nsim=-5))
try(spect(po,kernel.width=3,nsim=NA))
po@data[3] <- NA
try(spect(po,kernel.width=3,nsim=100))
po <- simulate(po,times=c(1:50,53:100))
try(spect(po,kernel.width=3,nsim=100))
po <- ricker
try(spect(pomp(po,rmeasure=function(x,t,params,...)
  c(y=if(x["N"]<10)rpois(n=1,lambda=x)else NA)),kernel.width=7,nsim=100))

try(spect.match(spp,weights="A"))
try(spect.match(spp,weights=rep(1,5)))
spp <- spect.match(spp,nsim=100,est="sigma",reltol=1e-3,maxit=100)
spp <- spect.match(spp,weights=function(freq)ifelse(freq<0.2,1,0),
                   reltol=1e-3,maxit=100,start=as.list(coef(spp)))
try(spect.match(spp,weights=function(freq)ifelse(freq<0.2,1,NA)))
try(spect.match(spp,est=spp@est,
                weights=function(freq)ifelse(freq<0.2,1,"C")))
try(spect.match(po,est="bob",nsim=100,kernel.width=3))
summary(spect.match(po,est="",nsim=10,kernel.width=3))$msg
summary(spect.match(po,est=NA,nsim=10,kernel.width=3))$msg
summary(spect.match(po,est=NULL,nsim=10,kernel.width=3))$msg
summary(spect.match(sp))$msg
try(summary(spect.match(po,est=NULL,vars="bob",nsim=10,kernel.width=3)))
try(summary(spect.match(po,est=NULL,nsim=NA,kernel.width=3)))
try(summary(spect.match(po,est=NULL,nsim=NULL,kernel.width=3)))
try(summary(spect.match(po,est=NULL,nsim=-3,kernel.width=3)))
try(summary(spect.match(po,est=NULL,kernel.width=3)))
try(summary(spect.match(po,est=NULL,nsim=20)))

dat <- as.data.frame(matrix(runif(60),20,3))
names(dat) <- letters[1:3]
dat$time <- 1:20
pomp(dat,times='time',t0=0,
     rprocess=euler.sim(Csnippet("
       x = rnorm(0,1);
       y = rnorm(0,1);
       z = rnorm(0,1);"),delta.t=1),
     rmeasure=Csnippet("
       a = rnorm(x,1);
       b = rnorm(y,1);
       c = rnorm(z,1);"),
     initializer=Csnippet("x = y = z = 0;"),
     statenames=tail(letters,3),
     params=c(dummy=1)) -> bob

plot(spect(bob,kernel.width=3,nsim=500),
     data.style=list(lwd=c(2,3),lty=2,col='red'))

## designed to fail
po <- ou2
coef(po,c("sigma.1","sigma.2","sigma.3","tau")) <- 0
spp <- spect.match(po,kernel.width=3,nsim=100,maxit=10,reltol=0.1,
                   est=c("alpha.1","alpha.2"))
stopifnot(is.na(spp@value))

dev.off()
