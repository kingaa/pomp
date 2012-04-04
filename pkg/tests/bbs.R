library(pomp)

data(bbs)

set.seed(48857734L)

coef(bbs)
coef(bbs,transform=TRUE)

prior <- parmat(coef(bbs),nrep=1000)
prior["log.beta1",] <- runif(n=1000,min=1,max=2)
prior["sigma",] <- runif(n=1000,min=2,max=4)
fit1 <- bsmc(bbs,params=prior,transform=TRUE,est=c("log.beta1","sigma"),smooth=0.2)
signif(coef(fit1),3)

fit2 <- traj.match(bbs,est=c("log.beta1","sigma"),transform=TRUE)
signif(coef(fit2),3)
