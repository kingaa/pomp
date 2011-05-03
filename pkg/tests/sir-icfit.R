library(pomp)

set.seed(343435488L)

pdf(file="sir-icfit.pdf")

data(euler.sir)
po <- window(euler.sir,end=0.25)
guess <- coef(po)
ics <- c("S.0","I.0","R.0") 
guess[ics[-3]] <- guess[ics[-3]]+c(0.5,-0.3)

plist <- list(
              probe.marginal("reports",ref=obs(po),order=3,diff=1,transform=sqrt),
              probe.acf("reports",lags=c(0,1,2,3,4,5),transform=sqrt),
              median=probe.median("reports")
              )

summary(pm.true <- probe(po,probes=plist,nsim=100,seed=1066L))

summary(pm.guess <- probe(po,params=guess,probes=plist,nsim=100,seed=1066L))

pm.fit <- probe.match(
                      po,
                      start=guess,
                      probes=plist,
                      est=ics[-1],
                      method="Nelder-Mead",
                      trace=3,
                      reltol=1e-5,
                      parscale=c(0.1,0.1),
                      nsim=100,
                      seed=1066L
                     )

summary(pm.fit)

comp.table <- cbind(true=exp(coef(po,ics)),guess=exp(guess[ics]),fit=exp(coef(pm.fit,ics)))
comp.table <- apply(comp.table,2,function(x)x/sum(x))
comp.table <- rbind(
                    comp.table,
                    synth.loglik=c(
                      summary(pm.true)$synth.loglik,
                      summary(pm.guess)$synth.loglik,
                      summary(pm.fit)$synth.loglik
                      )
                    )
comp.table

x <- sapply(
            list(true=pm.true,guess=pm.guess,fit=pm.fit),
            function (x) trajectory(x,times=time(x),t0=timezero(x))["cases",1,]
            )

plot(range(time(po)),range(c(states(po,"cases"),x)),bty='l',xlab="time",ylab="cases",type='n')
points(time(po),states(po,"cases"))
matlines(time(po),x,lty=1,col=c("red","blue","green"))
legend("topright",lty=1,bty='n',col=c("red","blue","green"),legend=colnames(x))

data(euler.sir)
po <- window(euler.sir,end=0.25)
guess <- coef(po)
ics <- c("S.0","I.0","R.0") 
guess[ics[-3]] <- guess[ics[-3]]+c(0.1,-0.2)

summary(tm.true <- traj.match(po,eval.only=TRUE))

summary(tm.guess <- traj.match(po,start=guess,eval.only=TRUE))

tm.fit <- traj.match(
                     po,
                     start=guess,
                     est=ics[-1],
                     method="sannbox",
                     maxit=300,
                     trace=2,
                     parscale=c(0.1,0.1)
                     )

tm.fit <- traj.match(
                     tm.fit,
                     est=ics[-1],
                     method="Nelder-Mead",
                     trace=3,
                     reltol=1e-8,
                     parscale=c(0.1,0.1)
                     )

summary(tm.fit)

comp.table <- cbind(true=exp(coef(po,ics)),guess=exp(guess[ics]),fit=exp(coef(tm.fit,ics)))
comp.table <- apply(comp.table,2,function(x)x/sum(x))
comp.table <- rbind(
                    comp.table,
                    loglik=sapply(list(tm.true,tm.guess,tm.fit),logLik)
                    )
comp.table

x <- sapply(
            list(true=tm.true,guess=tm.guess,fit=tm.fit),
            function (x) trajectory(x,times=time(x),t0=timezero(x))["cases",1,]
            )

plot(range(time(po)),range(c(states(po,"cases"),x)),bty='l',xlab="time",ylab="cases",type='n')
points(time(po),states(po,"cases"))
matlines(time(po),x,lty=1,col=c("red","blue","green"))
legend("topright",lty=c(NA,1,1,1),pch=c(1,NA,NA,NA),bty='n',col=c("black","red","blue","green"),legend=c("actual",colnames(x)))

dev.off()

