library(pomp.devel)

data(ou2)
true.p <- c(
            alpha.1=0.9,alpha.2=0,alpha.3=0,alpha.4=0.99,
            sigma.1=1,sigma.2=0,sigma.3=2,
            tau=1,x1.0=50,x2.0=-50
            )
simdata <- simulate(ou2,nsim=5,params=true.p,seed=394885)
guess.p <- true.p
guess.p[grep('sigma',names(guess.p))] <- 0

x <- sapply(
            simdata,
            function (d) {
              res <- traj.match(
                                d,
                                params=guess.p,
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
            bias=apply(x[names(true.p),],1,mean)-true.p,
            se=apply(x[names(true.p),],1,sd)
            )
      )
