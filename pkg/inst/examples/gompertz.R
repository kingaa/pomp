require(pomp)

## First code up the Gompertz example in R:

pomp(
     data=data.frame(time=1:100,Y=NA),
     times="time",
     t0=0,
     rprocess=discrete.time.sim( # a discrete-time process (see ?plugins)
       step.fun=function (x, t, params, delta.t, ...) { # this function takes one step t -> t+delta.t
         ## unpack and untransform the parameters:
         r <- exp(params["log.r"])
         K <- exp(params["log.K"])
         sigma <- exp(params["log.sigma"])
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
       ## unpack and untransform the parameters:
       tau <- exp(params["log.tau"])
       ## state at time t:
       X <- x["X"]
       ## generate a simulated observation:
       y <- c(Y=unname(rlnorm(n=1,meanlog=log(X),sdlog=tau)))
       return(y)
     },
     dmeasure=function (y, x, t, params, log, ...) { # measurement model density
       ## unpack and untransform the parameters:
       tau <- exp(params["log.tau"])
       ## state at time t:
       X <- x["X"]
       ## observation at time t:
       Y <- y["Y"]
       ## compute the likelihood of Y|X,tau
       f <- dlnorm(x=Y,meanlog=log(X),sdlog=tau,log=log)
       return(f)
     }
     ) -> gompertz


## Now code up the Gompertz example using native routines results in much faster computations.
## The C codes are included in the "examples" directory (file "gompertz.c")

po <- pomp(
           data=data.frame(time=seq(0,100,by=1),Y=NA),
           times="time",
           t0=0,
           rprocess=discrete.time.sim(
             step.fun="gompertz_simulator",
             PACKAGE="gompertz"
             ),
           rmeasure="gompertz_normal_rmeasure",
           dmeasure="gompertz_normal_dmeasure",
           skeleton.type="map",
           skeleton="gompertz_skeleton",
           PACKAGE="gompertz",
           paramnames=c("log.r","log.K","log.sigma","log.tau"),
           statenames=c("X"),
           obsnames=c("Y")
           )

params <- c(
            log.K=log(1),
            log.r=log(0.1),
            log.sigma=log(0.1),
            log.tau=log(0.1),
            X.0=1
            )

if (Sys.info()['sysname']=='Linux') {   # only run this under linux

  model <- "gompertz"
  pkg <- "pomp"
  modelfile <- paste(model,".c",sep="")
  solib <- paste(model,.Platform$dynlib.ext,sep="")

  ## compile the model into a shared-object library
  if (!file.copy(from=system.file(paste("examples/",modelfile,sep=""),package=pkg),to=getwd()))
    stop("cannot copy source code ",modelfile," to ",getwd())
  if (!file.copy(from=system.file("include/pomp.h",package=pkg),to=getwd()))
    stop("cannot copy header file ",modelfile," to ",getwd())
  cmd <- paste(R.home("bin/R"),"CMD SHLIB -o",solib,modelfile)
  rv <- system(cmd)
  if (rv!=0)
    stop("cannot compile shared-object library ",solib)

  dyn.load(solib) ## load the shared-object library

  ## compute a trajectory of the deterministic skeleton
  tic <- Sys.time()
  X <- trajectory(po,params=params)
  toc <- Sys.time()
  print(toc-tic)

  ## simulate from the model
  tic <- Sys.time()
  x <- simulate(po,params=params,nsim=3)
  toc <- Sys.time()
  print(toc-tic)
  
  dyn.unload(solib)

}
