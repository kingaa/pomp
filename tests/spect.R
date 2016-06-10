library(pomp)
pompExample(ou2)

set.seed(362083261L)

pdf(file="spect.pdf")

gm1 <- spect.match(ou2,
                  kernel.width=3,
                  detrend="mean",
                  nsim=50,
                  est=c("alpha.1","alpha.4"),
                  method="Nelder-Mead")
gm1@value
plot(gm1)

gm2 <- spect.match(ou2,
                   kernel.width=3,
                   detrend="mean",
                   nsim=49,
                   est=c("alpha.1","alpha.4"),
                   method="Nelder-Mead")
gm2@value
plot(gm2)

gm3 <- spect.match(ou2,
                   kernel.width=3,
                   detrend="linear",
                   nsim=50,
                   est=c("alpha.1","alpha.4"),
                   method="Nelder-Mead")
gm3@value
plot(gm3)

gm4 <- spect.match(ou2,
                   kernel.width=3,
                   detrend="quadratic",
                   nsim=50,
                   est=c("alpha.1","alpha.4"),
                   method="Nelder-Mead")
gm4@value
plot(gm4)

pompExample(ricker)

set.seed(6457673L)

sp <- spect(
            ricker,
            kernel.width=3,
            nsim=500,
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
            nsim=500,
            seed=838775L
            )
plot(sp)
invisible(summary(sp))

po <- ricker
coef(po,"phi") <- 30
sp <- spect(
            po,
            kernel.width=3,
            nsim=500,
            seed=838775L
            )
plot(sp)
invisible(summary(sp))

plot(simulate(sp),variables="y")

dev.off()
