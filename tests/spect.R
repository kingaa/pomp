library(pomp)
pompExample(ou2)

set.seed(362083261L)

png(filename="spect-%02d.png",res=100)

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
                   method="Nelder-Mead")
plot(gm4)

pompExample(ricker)

sp <- spect(ricker,kernel.width=3,nsim=100,seed=838775L)
invisible(summary(sp))

spp <- spect.match(sp,eval.only=TRUE)
invisible(summary(spp))

spp <- spect.match(sp,nsim=100,est=c("sigma","phi"),reltol=1e-3,maxit=100)
invisible(summary(spp))

po <- ricker
coef(po,"r") <- 5
sp <- spect(po,kernel.width=1,nsim=100,seed=838775L)
invisible(summary(sp))

po <- ricker
coef(po,"phi") <- 30
try(spect(po,nsim=100,seed=838775L))
sp <- spect(po,kernel.width=11,nsim=100,seed=838775L)
invisible(summary(sp))

plot(simulate(sp),variables="y")

try(spect.match(spp,est=spp@est,weights="A"))
try(spect.match(spp,est=spp@est,weights=rep(1,5)))
spp <- spect.match(spp,est=spp@est,
                   weights=function(freq)ifelse(freq<0.2,1,0))
try(spect.match(spp,est=spp@est,
                weights=function(freq)ifelse(freq<0.2,1,NA)))
try(spect.match(spp,est=spp@est,
                weights=function(freq)ifelse(freq<0.2,1,"C")))


dev.off()
