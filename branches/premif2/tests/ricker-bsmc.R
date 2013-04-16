library(pomp)

data(ricker)

pdf(file="ricker-bsmc.pdf")

set.seed(6457673L)

po <- ricker

Np <- 10000
params <- parmat(coef(ricker),nrep=Np)
params["r",] <- exp(runif(n=Np,min=2,max=5))
params["sigma",] <- runif(n=Np,min=0.1,max=1)

fit <- bsmc(ricker,params=params,est=c("r","sigma"),transform=TRUE,smooth=0.2)

invisible(apply(fit$prior[c("r","sigma"),],1,mean))

invisible(apply(fit$post[c("r","sigma"),],1,mean))

invisible(coef(fit))

plot(fit,thin=300)

dev.off()
