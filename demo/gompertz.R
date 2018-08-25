library(pomp)
library(ggplot2)
library(magrittr)

## First, code up the Gompertz example in R:

simulate(times=1:100,t0=0,
  params=c(K=1,r=0.1,sigma=0.1,tau=0.1,X.0=1),
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
  dmeasure=function (tau, X, Y, ..., log) {
    dlnorm(x=Y,meanlog=log(X),sdlog=tau,log=log)
  },
  partrans=parameter_trans(
    toEst=function(params,...){
      log(params)
    },
    fromEst=function(params,...){
      exp(params)
    }
  )
) -> gompertz

## Now code up the Gompertz example using C snippets: results in much faster computations.

simulate(times=0:100,t0=0,
  params=c(K=1,r=0.1,sigma=0.1,tau=0.1,X.0=1),
  dmeasure=Csnippet("
    lik = dlnorm(Y,log(X),tau,give_log);"
  ),
  rmeasure=Csnippet("
    Y = rlnorm(log(X),tau);"
  ),
  rprocess=discrete.time.sim(
    step.fun=Csnippet("
    double S = exp(-r*dt);
    double logeps = (sigma > 0.0) ? rnorm(0,sigma) : 0.0;
    /* note that X is over-written by the next line */
    X = pow(K,(1-S))*pow(X,S)*exp(logeps);"
    ),
    delta.t=1
  ),
  skeleton=map(Csnippet("
    double dt = 1.0;
    double S = exp(-r*dt);
    /* note that X is not over-written in the skeleton function */
    DX = pow(K,1-S)*pow(X,S);"
  ),delta.t=1),
  partrans=parameter_trans(log=c("r","K","sigma","tau","X.0")),
  paramnames=c("r","K","sigma","X.0","tau"),
  obsnames="Y",
  statenames="X"
) -> Gompertz

p <- parmat(coef(Gompertz),nrep=4)
p["X.0",] <- c(0.5,0.9,1.1,1.5)

## compute a trajectory of the deterministic skeleton
X <- trajectory(Gompertz,params=p,format="data.frame")
X %>%
  ggplot(aes(x=time,y=X,group=.id,color=.id))+
  guides(color=FALSE)+
  geom_line()+
  theme_bw()+
  labs(title="Gompertz model",subtitle="deterministic trajectories")

## simulate from the model
x <- simulate(Gompertz,params=p,format="data.frame")

x %>%
  ggplot(aes(x=time,y=X,group=.id,color=.id))+
  geom_line()+
  guides(color=FALSE)+
  theme_bw()+
  labs(title="Gompertz model",subtitle="stochastic simulations")

## run a particle filter
pf <- replicate(n=10,pfilter(Gompertz,Np=500))

print(logmeanexp(sapply(pf,logLik),se=TRUE))
