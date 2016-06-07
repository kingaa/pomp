library(pomp)

pdf(file="ricker-spect.pdf")

pompExample(ricker)

set.seed(6457673L)

sp <- spect(
            ricker,
            kernel.width=3,
            nsim=1000,
            seed=838775L
            )
plot(sp)
invisible(summary(sp))

spp <- spect.match(sp,eval.only=TRUE)
plot(spp)
invisible(summary(spp))

spp <- spect.match(sp,nsim=100,est=c("sigma","phi"))
plot(spp)
invisible(summary(spp))

po <- ricker
coef(po,"r") <- 5
sp <- spect(
            po,
            kernel.width=3,
            nsim=1000,
            seed=838775L
            )
plot(sp)
invisible(summary(sp))

po <- ricker
coef(po,"phi") <- 30
sp <- spect(
            po,
            kernel.width=3,
            nsim=1000,
            seed=838775L
            )
plot(sp)
invisible(summary(sp))

plot(simulate(sp),variables="y")

dev.off()
