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
     parameter.inv.transform=function(params,...){
       params <- log(params[c("X.0","r","K","tau","sigma")])
       names(params) <- c("log.X.0","log.r","log.K","log.tau","log.sigma")
       params
     },
     parameter.transform=function(params,...){
       params <- exp(params[c("log.X.0","log.r","log.K","log.tau","log.sigma")])
       names(params) <- c("X.0","r","K","tau","sigma")
       params
     }
     ) -> gompertz

## Now code up the Gompertz example using native routines results in much faster computations.
## The C codes are included in the "examples" directory (file "gompertz.c")

if (Sys.info()['sysname']=='Linux') {   # only run this under linux

  model <- "gompertz"
  ## name of the file holding the native codes for this model:
  modelfile <- system.file("examples",paste(model,".c",sep=""),package="pomp")
  ## name of the shared-object library
  solib <- paste(model,.Platform$dynlib.ext,sep="")
  ## environment variables needed to locate the pomp header file:
  cflags <- paste("PKG_CFLAGS=\"-I",system.file("include/",package="pomp"),"\"")

  ## compile the shared-object library containing the model codes:
  rv <- system2(
                command=R.home("bin/R"),
                args=c("CMD","SHLIB","-o",solib,modelfile),
                env=cflags
                )
  ## compile the shared-object library containing the model codes:
  if (rv!=0)
    stop("cannot compile shared-object library ",sQuote(solib))

  dyn.load(solib) ## load the shared-object library

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
             paramnames=c("r","K","sigma","tau"),
             statenames=c("X"),
             obsnames=c("Y"),
             parameter.transform=function(params,...){
               params <- log(params[c("X.0","r","K","tau","sigma")])
               names(params) <- c("log.X.0","log.r","log.K","log.tau","log.sigma")
               params
             },
             parameter.inv.transform=function(params,...){
               params <- exp(params[c("log.X.0","log.r","log.K","log.tau","log.sigma")])
               names(params) <- c("X.0","r","K","tau","sigma")
               params
             }
             )

  ## set the parameters
  coef(po) <- c(K=1,r=0.1,sigma=0.1,tau=0.1,X.0=1)

  p <- parmat(coef(po),nrep=4)
  p["X.0",] <- c(0.5,0.9,1.1,1.5)
  
  ## compute a trajectory of the deterministic skeleton
  tic <- Sys.time()
  X <- trajectory(po,params=p,as.data.frame=TRUE)
  toc <- Sys.time()
  print(toc-tic)
  X <- reshape(X,dir="wide",v.names="X",timevar="traj",idvar="time")
  matplot(X$time,X[-1],type='l',lty=1,bty='l',xlab="time",ylab="X",
          main="Gompertz model\ndeterministic trajectories")
  
  ## simulate from the model
  tic <- Sys.time()
  x <- simulate(po,params=p,as.data.frame=TRUE)
  toc <- Sys.time()
  print(toc-tic)

  x <- reshape(x,dir="wide",v.names=c("Y","X"),timevar="sim",idvar="time")
  op <- par(mfrow=c(2,1),mgp=c(2,1,0),mar=c(3,3,0,0),bty='l')
  matplot(x$time,x[c("X.1","X.2","X.3")],lty=1,type='l',xlab="time",ylab="X",
          main="Gompertz model\nstochastic simulations")
  matplot(x$time,x[c("Y.1","Y.2","Y.3")],lty=1,type='l',xlab="time",ylab="Y")
  par(op)
  
  dyn.unload(solib)

}
