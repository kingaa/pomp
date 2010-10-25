library(pomp)

data(ricker)

pdf(file="ricker.pdf")

tj.1 <- trajectory(ricker)
plot(time(ricker),tj.1[1,,],type='l')
tj.2 <- trajectory(ricker,times=c(30:50),t0=0)
lines(30:50,tj.2[1,,],col='red',lwd=2)
max(abs(tj.1[,,time(ricker)>=30]-tj.2[,,]))

po <- ricker
try(
    coef(po,"r")
    )
coef(po,c("r","phi")) <- c(0,0)
coef(po,c("log.r","log.phi")) <- c(a=0,b=0)
coef(po,c("log.r","log.phi")) <- 0
coef(po) <- c(log.phi=0,log.r=3.5,N.0=10,e.0=0,log.sigma=-Inf)
coef(po)
coef(po,"new") <- 3
plot(simulate(po))
coef(po)

set.seed(64857673L)

guess <- ricker
coef(guess) <- c(log.r=log(20),log.sigma=log(1),log.phi=log(20),N.0=7,e.0=0)

mf <- mif(
          guess,
          Nmif=100,
          Np=1000,
          cooling.factor=0.99,
          var.factor=2,
          ic.lag=3,
          max.fail=50,
          rw.sd=c(log.r=0.1,log.sigma=0.1,log.phi=0.1)
          )
mf <- continue(mf,Nmif=500,max.fail=20)

plist <- list(
              probe.marginal("y",ref=obs(ricker),transform=sqrt),
              probe.acf("y",lags=c(0,1,2,3,4),transform=sqrt),
              probe.nlar("y",lags=c(1,1,1,2),powers=c(1,2,3,1),transform=sqrt)
              )

pr.guess <- probe(guess,probes=plist,nsim=1000,seed=1066L)
pr.tr <- probe(ricker,probes=plist,nsim=1000,seed=1066L)

pm <- probe.match(
                  pr.guess,
                  est=c("log.r","log.sigma","log.phi"),
                  method="Nelder-Mead",
                  maxit=2000,
                  seed=1066L,
                  reltol=1e-8,
                  trace=3
                  )

pf.tr <- pfilter(ricker,Np=1000,max.fail=50,seed=1066L)
pf.guess <- pfilter(guess,Np=1000,max.fail=50,seed=1066L)
pf.mf <- pfilter(mf,Np=1000,seed=1066L)
pf.pm <- pfilter(pm,Np=1000,max.fail=10,seed=1066L)

pr.mf <- probe(mf,nsim=1000,probes=plist,seed=1066L)

res <- rbind(
             cbind(guess=coef(guess),truth=coef(ricker),MLE=coef(mf),PM=coef(pm)),
             loglik=c(
               pf.guess$loglik,
               pf.tr$loglik,
               pf.mf$loglik,
               pf.pm$loglik
               ),
             synth.loglik=c(
               summary(pr.guess)$synth.loglik,
               summary(pr.tr)$synth.loglik,
               summary(pr.mf)$synth.loglik,
               summary(pm)$synth.loglik
               )
             )

print(res,digits=3)

plot(ricker)
plot(simulate(guess))
plot(simulate(mf))
plot(simulate(pm))

dev.off()
