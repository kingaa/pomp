require(pomp)

## a simple two-dimensional random walk
## this makes use of the 'onestep.sim' plugin
## which we can use since we can simulate the
## increment of a random walk over any time

rw2 <- pomp(
            rprocess = onestep.sim(
              step.fun = function(x, t, params, delta.t, ...) {
                c(
                  x1=rnorm(n=1,mean=x['x1'],sd=params['s1']*delta.t),
                  x2=rnorm(n=1,mean=x['x2'],sd=params['s2']*delta.t)
                  )
              }
              ),
            dprocess = onestep.dens(
              dens.fun = function (x1, t1, x2, t2, params, ...) {
                sum(
                    dnorm(
                          x=x2[c('x1','x2')],
                          mean=x1[c('x1','x2')],
                          sd=params[c('s1','s2')]*(t2-t1),
                          log=TRUE
                          ),
                    na.rm=TRUE
                    )
              }
              ),
            measurement.model=list(
              y1 ~ norm(mean=x1,sd=tau),
              y2 ~ norm(mean=x2,sd=tau)
              ),
            times=1:100,
            data=rbind(
              y1=rep(0,100),
              y2=rep(0,100)
              ),
            t0=0
            )

rw2 <- simulate(rw2,params=c(s1=3,s2=1,x1.0=0,x2.0=1,tau=10))
plot(rw2)
