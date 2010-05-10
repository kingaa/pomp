library(pomp)

data(ou2)
true.p <- coef(ou2)
simdata <- simulate(ou2,nsim=5,params=true.p,seed=394885)
guess.p <- true.p
guess.p[grep('sigma',names(guess.p))] <- 0

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
              c(conv=res$convergence,loglik=res$value,res$params)
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
