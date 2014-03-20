library(pomp)

set.seed(93885485L)

pompExample(ou2)
true.p <- coef(ou2)
simdata <- simulate(ou2,nsim=5,params=true.p,seed=394885)
guess.p <- true.p

x <- sapply(
            simdata,
            function (d) {
              res <- traj.match(
                                d,
                                start=guess.p,
                                est=c('alpha.1','alpha.4','x1.0','x2.0','tau'),
                                maxit=2000,
                                reltol=1e-8
                                )
              c(conv=res$convergence,loglik=logLik(res),coef(res))
            }
            )
range(x['conv',])
range(x['loglik',])
print(
      cbind(
            truth=true.p,
            mean.est=apply(x[names(true.p),],1,mean),
            bias=apply(x[names(true.p),],1,mean)-true.p,
            se=apply(x[names(true.p),],1,sd)
            ),
      digits=4
      )

summary(traj.match(ou2,est=c('alpha.1','alpha.4','x1.0','x2.0','tau'),method="subplex",maxit=100))

summary(traj.match(ou2,est=c('alpha.1','x1.0','alpha.4','x2.0','tau'),method="sannbox",trace=1,parscale=0.1,maxit=100))

summary(traj.match(ou2))

ofun <- traj.match.objfun(ou2,est=c('x1.0','x2.0','alpha.1','alpha.4','tau'))
print(ofun(coef(ou2,c('x1.0','x2.0','alpha.1','alpha.4','tau'))))

try(ofun <- traj.match.objfun(ou2,est=c('x1.0','x2.0','bob','alpha.4','harry')))
try(traj.match(ou2,est=c('x1.0','x2.0','bob','alpha.4','tau')))

fit <- optim(par=c(0,0,1,0.5,0.5),fn=ofun,method="Nelder-Mead",control=list(maxit=400))
print(fit)
stopifnot(fit$convergence==0)

pompExample(ou2)
p <- coef(ou2)
ou2@params <- numeric(0)
res <- traj.match(ou2,start=p)
print(coef(res),digits=4)
res <- traj.match(res,est=c('alpha.1','alpha.4','x1.0','x2.0','tau'),maxit=2000,reltol=1e-8)
print(coef(res),digits=4)
print(p,digits=4)
