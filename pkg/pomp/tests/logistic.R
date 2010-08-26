library(pomp)

po <- pomp(
           data=rbind(obs=rep(0,1000)),
           times=0.1*seq.int(length=1000),
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
                          sd=sigma*n*sqrt(delta.t)
                          )
                    )
             }
             ),
           measurement.model=obs~lnorm(meanlog=log(n),sdlog=log(1+tau)),
           skeleton.vectorfield=function(x,t,params,...){
             with(
                  as.list(c(x,params)),
                  r*n*(1-n/K)
                  )
           }
           )

params <- c(n.0=10000,K=10000,r=0.9,sigma=0.4,tau=0.1)
set.seed(73658676)
po <- simulate(po,params=params)

t <- seq(0,by=0.005,length=50)
x <- simulate(po,times=t,states=T,obs=T)

print(
      dprocess(
               po,
               x=x$states[,,45:50,drop=F],
               times=t[45:50],
               params=as.matrix(params),
               log=TRUE
               ),
      digits=4
      )

print(
      dmeasure(
               po,
               y=rbind(obs=x$obs[,1,45:50]),
               x=x$states[,,45:50,drop=F],
               times=t[45:50],
               params=as.matrix(params),
               log=TRUE
               ),
      digits=4
      )

print(
      drop(
           skeleton(
                    po,
                    x=array(
                      seq(0,12000,by=1000),
                      dim=c(1,1,13),
                      dimnames=list('n',NULL,NULL)
                      ),
                    t=rep(0,13),
                    params=as.matrix(params)
                    )
           ),
      digits=4
      )
               
params <- cbind(c(n.0=100,K=10000,r=0.2,sigma=0.4,tau=0.1),c(n.0=1000,K=11000,r=0.1,sigma=0.4,tau=0.1))
x <- trajectory(po,params=params)
pdf(file='logistic.pdf')
plot(po)
matplot(time(po,t0=TRUE),t(x['n',,]),type='l')
dev.off()
