library(pomp)

data(ricker)

x <- array(
           data=states(ricker),
           dim=c(2,3,52),
           dimnames=list(rownames(states(ricker)),NULL,NULL)
           )
p <- array(
           data=coef(ricker),
           dim=c(5,3),
           dimnames=list(names(coef(ricker)),NULL)
           )
p["log.r",] <- c(1,2,4)
f <- skeleton(ricker,x=x,params=p,t=time(ricker,t0=T))
plot(x[1,,],f[1,,],type='n')
points(x[1,1,],f[1,1,],col='red')
points(x[1,2,],f[1,2,],col='blue')
points(x[1,3,],f[1,3,],col='green')

## non-autonomous case

data(euler.sir)
x <- array(
           data=states(euler.sir)[,1],
           dim=c(5,3,100),
           dimnames=list(rownames(states(euler.sir)),NULL,NULL)
           )
p <- array(
           data=coef(euler.sir),
           dim=c(15,3),
           dimnames=list(names(coef(euler.sir)),NULL)
           )
p["log.beta2",1:2] <- c(3,5)  ## try different values of one of the seasonality parameters
## compute the skeleton at each point
f <- skeleton(euler.sir,x=x,params=p,t=seq(0,1,length=100))
## verify that the skeleton varies with time
matplot(seq(0,1,length=100),t(f[1,,]),type='l',lty=1)

x <- trajectory(euler.sir) ## use deSolve to compute a deterministic trajectory
f <- skeleton(euler.sir,x=x,params=p[,3,drop=F],t=time(euler.sir))  ## evaluate skeleton at each point along it

fit <- smooth.spline(time(euler.sir),x["S",,])## fit a spline to the trajectory
plot(predict(fit,deriv=1)$y,f["S",,])## compare spline estimate with computed vectorfield
abline(a=0,b=1)

fit <- smooth.spline(time(euler.sir),x["I",,])
plot(predict(fit,deriv=1)$y,f["I",,])
abline(a=0,b=1)

pomp.skeleton <- function(times,y,p,more) {
# Turns a skeleton function from a 'pomp' object into the right hand
# side of and ODE for use in CollocInfer

  x <- array(y,dim=c(nrow(y),1,ncol(y)),dimnames=list(rownames(y),NULL,NULL))
  params <- array(data=p,dim=c(length(p),1),dimnames=list(names(p),NULL))

  skeleton(more$pomp.obj,x=x,params=params,t=times)
}

x <- array(
           data=states(ricker),
           dim=c(2,3,51),
           dimnames=list(rownames(states(ricker)),NULL,NULL)
           )
p <- array(
           data=coef(ricker),
           dim=c(5,3),
           dimnames=list(names(coef(ricker)),NULL)
           )
p["log.r",]<- c(1,2,4)

f <- skeleton(ricker,x=x,params=p,t=time(ricker))

pomp.skeleton(time(ricker),x[,1,],p[,1],list(pomp.obj=ricker))

