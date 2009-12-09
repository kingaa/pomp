library(pomp)

set.seed(921625222L)

pf <- pfilter(ou2,Np=1000,save.states=TRUE)
ll <- cumsum(pf$cond.loglik)
pp <- matrix(data=coef(ou2),nrow=length(coef(ou2)),ncol=1000,dimnames=list(names(coef(ou2)),NULL))
for (t in seq(60,90,by=10)) {
  pp[c("x1.0","x2.0"),] <- pf$states[c("x1","x2"),,t]
  y <- simulate(ou2,params=pp,obs=TRUE,times=time(ou2)[t:(t+10)])
  mn <- apply(y,c(1,3),mean)
  sd <- apply(y,c(1,3),sd)
  z <- (data.array(ou2)[,t:(t+10)]-mn)/sd       ## z score
  mse <- (data.array(ou2)[,t:(t+10)]-mn)^2+sd^2 ## mean squared error
}
fit <- mif(ou2,Nmif=3,rw.sd=c(alpha.1=0.1,alpha.4=0.1),Np=1000,cooling.factor=0.98,var.factor=1,ic.lag=2)
try(
    pfilter(fit,save.states=TRUE)
    )

