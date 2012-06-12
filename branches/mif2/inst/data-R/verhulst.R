require(pomp)

simulate(
         pomp(
              data=rbind(obs=rep(0,1000)),
              times=seq(0.1,by=0.1,length=1000),
              t0=0,
              rprocess=euler.sim(
                step.fun=function(x,t,params,delta.t,...){
                  r <- params["r"]
                  K <- params["K"]
                  sigma <- params["sigma"]
                  n <- x["n"]
                  rnorm(
                        n=1,
                        mean=n+r*n*(1-n/K)*delta.t,
                        sd=sigma*n*sqrt(delta.t)
                        )
                },
                delta.t=0.01
                ),
              dprocess=onestep.dens(
                dens.fun=function(x1,x2,t1,t2,params,log,...){
                  delta.t <- t2-t1
                  r <- params["r"]
                  K <- params["K"]
                  sigma <- params["sigma"]
                  n <- x1["n"]
                  dnorm(
                        x=x2["n"],
                        mean=n+r*n*(1-n/K)*delta.t,
                        sd=sigma*n*sqrt(delta.t),
                        log=log
                        )
                }
                ),
              measurement.model=obs~lnorm(meanlog=log(n),sdlog=log(1+tau)),
              skeleton.type="vectorfield",
              skeleton=function(x,t,params,...){
                r <- params["r"]
                K <- params["K"]
                n <- x["n"]
                f <- r*n*(1-n/K)
                names(f) <- "n"
                f
              }
              ),
         params=c(
           n.0=10000,
           K=10000,
           r=0.9,
           sigma=0.4,
           tau=0.1
           ),
         nsim=1,
         seed=73658676L
         ) -> verhulst

save(verhulst,file="verhulst.rda",compress="xz")
