require(pomp)

## First, code up the Gompertz example in R:

pomp(
     data=data.frame(time=1:100,Y=NA),
     times="time",
     t0=0,
     rprocess=discrete.time.sim( # a discrete-time process (see ?plugins)
       step.fun=function (x, t, params, delta.t, ...) { # this function takes one step t -> t+delta.t
         ## unpack the parameters:
         r <- params["r"]
         K <- params["K"]
         sigma <- params["sigma"]
         ## the state at time t:
         X <- x["X"]
         ## generate a log-normal random variable:
         eps <- exp(rnorm(n=1,mean=0,sd=sigma))
         ## compute the state at time t+delta.t:
         S <- exp(-r*delta.t)
         xnew <- c(X=unname(K^(1-S)*X^S*eps))
         return(xnew)
       },
       delta.t=1                  # the size of the discrete time-step
       ),
     rmeasure=function (x, t, params, ...) {# the measurement model simulator
       ## unpack the parameters:
       tau <- params["tau"]
       ## state at time t:
       X <- x["X"]
       ## generate a simulated observation:
       y <- c(Y=unname(rlnorm(n=1,meanlog=log(X),sdlog=tau)))
       return(y)
     },
     dmeasure=function (y, x, t, params, log, ...) { # measurement model density
       ## unpack the parameters:
       tau <- params["tau"]
       ## state at time t:
       X <- x["X"]
       ## observation at time t:
       Y <- y["Y"]
       ## compute the likelihood of Y|X,tau
       f <- dlnorm(x=Y,meanlog=log(X),sdlog=tau,log=log)
       return(f)
     },
     toEstimationScale=function(params,...){
       log(params)
     },
     fromEstimationScale=function(params,...){
       exp(params)
     }
     ) -> gompertz

## Now code up the Gompertz example using C snippets: results in much faster computations.

dmeas <- "
    lik = dlnorm(Y,log(X),tau,give_log);
"

rmeas <- "
    Y = rlnorm(log(X),tau);
"

step.fun <- "
  double S = exp(-r*dt);
  double logeps = (sigma > 0.0) ? rnorm(0,sigma) : 0.0;
  /* note that X is over-written by the next line */
  X = pow(K,(1-S))*pow(X,S)*exp(logeps); 
"

skel <- "
  double dt = 1.0;
  double S = exp(-r*dt);
  /* note that X is not over-written in the skeleton function */
  DX = pow(K,1-S)*pow(X,S); 
"

partrans <- "
  Tr = log(r);
  TK = log(K);
  Tsigma = log(sigma);
  TX_0 = log(X_0);
  Ttau = log(tau);
"

paruntrans <- "
  Tr = exp(r);
  TK = exp(K);
  Tsigma = exp(sigma);
  TX_0 = exp(X_0);
  Ttau = exp(tau);
"

pomp(
     data=data.frame(t=1:100,Y=NA),
     times="t",
     t0=0,
     paramnames=c("r","K","sigma","X.0","tau"),
     statenames=c("X"),
     dmeasure=Csnippet(dmeas),
     rmeasure=Csnippet(rmeas),
     rprocess=discrete.time.sim(
       step.fun=Csnippet(step.fun),
       delta.t=1
       ),
     skeleton=Csnippet(skel),
     skeleton.type="map",
     skelmap.delta.t=1,
     toEstimationScale=Csnippet(partrans),
     fromEstimationScale=Csnippet(paruntrans)
     ) -> Gompertz

## simulate some data
Gompertz <- simulate(Gompertz,params=c(K=1,r=0.1,sigma=0.1,tau=0.1,X.0=1))

p <- parmat(coef(Gompertz),nrep=4)
p["X.0",] <- c(0.5,0.9,1.1,1.5)

## compute a trajectory of the deterministic skeleton
tic <- Sys.time()
X <- trajectory(Gompertz,params=p,as.data.frame=TRUE)
toc <- Sys.time()
print(toc-tic)
X <- reshape(X,dir="wide",v.names="X",timevar="traj",idvar="time")
matplot(X$time,X[-1],type='l',lty=1,bty='l',xlab="time",ylab="X",
        main="Gompertz model\ndeterministic trajectories")

## simulate from the model
tic <- Sys.time()
x <- simulate(Gompertz,params=p,as.data.frame=TRUE)
toc <- Sys.time()
print(toc-tic)

x <- reshape(x,dir="wide",v.names=c("Y","X"),timevar="sim",idvar="time")
op <- par(mfrow=c(2,1),mgp=c(2,1,0),mar=c(3,3,0,0),bty='l')
matplot(x$time,x[c("X.1","X.2","X.3")],lty=1,type='l',xlab="time",ylab="X",
        main="Gompertz model\nstochastic simulations")
matplot(x$time,x[c("Y.1","Y.2","Y.3")],lty=1,type='l',xlab="time",ylab="Y")
par(op)

## run a particle filter
tic <- Sys.time()
pf <- replicate(n=10,pfilter(Gompertz,Np=500))
toc <- Sys.time()
print(toc-tic)

print(logmeanexp(sapply(pf,logLik),se=TRUE))
