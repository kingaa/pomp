library(pomp)

data(ricker)

pdf(file="ricker-spect.pdf")

set.seed(6457673L)

sp <- spect(
            ricker,
            kernel.width=3,
            nsim=1000,
            seed=838775L
            )
plot(sp)
summary(sp)

po <- ricker
coef(po,"log.r") <- log(5)
sp <- spect(
            po,
            kernel.width=3,
            nsim=1000,
            seed=838775L
            )
plot(sp)
summary(sp)

po <- ricker
coef(po,"log.phi") <- log(30)
sp <- spect(
            po,
            kernel.width=3,
            nsim=1000,
            seed=838775L
            )
plot(sp)
summary(sp)

dev.off()
