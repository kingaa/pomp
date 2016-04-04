library(pomp)

pompExample(bbs)

guess <- c(
           mu=0,gamma=1/3,beta=1,beta.sd=0,iota=0,
           pop=1400,rho=0.9,sigma=3.6,
           S_0=1390,I_0=1,R_0=0
           )
est <- c("beta","gamma")

tj1 <- trajectory(bbs,params=guess,as.data.frame=TRUE)
tail(tj1)

tj2 <- trajectory(bbs,params=guess,hmax=0.001,as.data.frame=TRUE)
tail(tj2)

tm1 <- traj.match(bbs,start=guess,transform=TRUE,est=est,method="subplex",reltol=1e-7)

tmf <- traj.match.objfun(bbs,params=guess,est=est,transform=TRUE,hmax=0.001)

library(subplex)
fit <- subplex(fn=tmf,par=log(guess[est]),control=list(reltol=1e-7))
tm2 <- bbs
coef(tm2) <- guess
coef(tm2,names(fit$par),transform=T) <- fit$par

round(coef(tm1,est)/coef(tm2,est),5)

options(warn=2)
ofun <- traj.match.objfun(window(bbs,end=3),est=c("beta","gamma"),transform=TRUE,maxsteps=10,rtol=1e-6)
try(optim(fn=ofun,par=c(0,-1),method="Nelder-Mead",control=list(reltol=1e-10)))
