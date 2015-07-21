library(pomp)

pompExample(ricker)

pdf(file="ricker-bsmc.pdf")

set.seed(6457673L)

po <- pomp(
           ricker,
           rprior=function (params, ...) {
             params["r"] <- exp(runif(n=1,min=2,max=5))
             params["sigma"] <- runif(n=1,min=0.1,max=1)
             params
           }
           )

Np <- 10000

fit <- bsmc2(po,Np=1000,est=c("r","sigma"),transform=TRUE,smooth=0.2)

invisible(apply(fit$prior[c("r","sigma"),],1,mean))

invisible(apply(fit$post[c("r","sigma"),],1,mean))

invisible(coef(fit))

plot(fit,thin=300)

dev.off()
