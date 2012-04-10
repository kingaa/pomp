require(pomp)

## a stochastic version of the Verhulst-Pearl logistic model
## this evolves in continuous time, but we approximate the
## stochastic dynamics using an Euler approximation
## (plugin 'euler.sim') with fixed step-size

po <- pomp(
           data=rbind(obs=rep(0,1000)),
           times=seq(0.1,by=0.1,length=1000),
           t0=0,
           rprocess=euler.sim(
             step.fun=function(x,t,params,delta.t,...){
               with(
                    as.list(c(x,params)),
                    rnorm(
                          n=1,
                          mean=n+r*n*(1-n/K)*delta.t,
                          sd=sigma*n*sqrt(delta.t)
                          )
                    )
             },
             delta.t=0.01
             ),
           dprocess=onestep.dens(
             dens.fun=function(x1,x2,t1,t2,params,log,...){
               delta.t <- t2-t1
               with(
                    as.list(c(x1,params)),
                    dnorm(
                          x=x2['n'],
                          mean=n+r*n*(1-n/K)*delta.t,
                          sd=sigma*n*sqrt(delta.t),
                          log=log
                          )
                    )
             }
             ),
           measurement.model=obs~lnorm(meanlog=log(n),sdlog=log(1+tau)),
           skeleton.type="vectorfield",
           skeleton=function(x,t,params,...){
             with(
                  as.list(c(x,params)),
                  r*n*(1-n/K)
                  )
           }
           )

params <- c(n.0=10000,K=10000,r=0.9,sigma=0.4,tau=0.1)
set.seed(73658676)
po <- simulate(po,params=params)
plot(po)

params <- cbind(
                c(n.0=100,K=10000,r=0.2,sigma=0.4,tau=0.1),
                c(n.0=1000,K=11000,r=0.1,sigma=0.4,tau=0.1)
                )
x <- trajectory(po,params=params,as.data.frame=TRUE)
x <- reshape(x,dir="wide",idvar="time",timevar="traj")
matplot(x$time,x[-1],type='l',bty='l',lty=1,xlab="time",ylab="n")
