library(pomp)

data(ricker)

pdf(file="ricker-bsmc.pdf")

set.seed(6457673L)

po <- ricker

Np <- 10000
params <- parmat(coef(ricker),nrep=Np)
params["log.r",] <- runif(n=Np,min=2,max=5)
params["sigma",] <- runif(n=Np,min=0.1,max=1)

fit <- bsmc(ricker,params=params,est=c("log.r","log.sigma"),transform=TRUE,smooth=0.2)

print(apply(fit$prior[c("log.r","log.sigma"),],1,mean))

print(apply(fit$post[c("log.r","log.sigma"),],1,mean))

plot(fit,breaks=30,thin=300)

dev.off()
