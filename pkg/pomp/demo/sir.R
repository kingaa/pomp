require(pomp)

## Coding up the SIR example using native routines results in much faster computations.
## The C codes "sir_euler_simulator", "sir_ODE", "binom_rmeasure", and "binom_dmeasure"
## are included in the "examples" directory (file "sir.c")

if (Sys.info()['sysname']=='Linux') {   # only run this under linux

  model <- "sir"
  ## name of the file holding the model codes:
  modelfile <- system.file("examples",paste(model,".c",sep=""),package="pomp")
  ## name of the shared-object library:
  solib <- paste(model,.Platform$dynlib.ext,sep="")
  ## environment variables needed to locate the pomp header file:
  cflags <- paste("PKG_CFLAGS=\"-I",system.file("include/",package="pomp"),"\"")

  ## compile the shared-object library containing the model codes:
  rv <- system2(
                command=R.home("bin/R"),
                args=c("CMD","SHLIB","-o",solib,modelfile),
                env=cflags
                )
  if (rv!=0)
    stop("cannot compile shared-object library ",sQuote(solib))

  po <- pomp(
             data=data.frame(
               time=seq(from=1/52,to=4,by=1/52),
               reports=NA
               ),
             times="time",
             t0=0,
             rprocess=euler.sim(
               step.fun="sir_euler_simulator", # native routine for the simulation step
               delta.t=1/52/20                 # Euler stepsize: 1/20 wk
               ),
             skeleton.type="vectorfield",
             skeleton="sir_ODE", # native routine for the skeleton
             rmeasure="binomial_rmeasure", # binomial measurement model
             dmeasure="binomial_dmeasure", # binomial measurement model
             PACKAGE="sir", ## name of the shared-object library
             ## the order of the observables assumed in the native routines:
             obsnames="reports",
             ## the order of the state variables assumed in the native routines:
             statenames=c("S","I","R","cases","W"),
             ## the order of the parameters assumed in the native routines:
             paramnames=c(
               "gamma","mu","iota","beta1","beta.sd","pop","rho",
               "nbasis","degree","period"
               ), 
             ## reset cases to zero after each new observation:
             zeronames=c("cases"),      
             logvar=c(
               "gamma","mu","iota",
               "beta1","beta2","beta3","beta.sd",
               "S.0","I.0","R.0"
               ),
             logitvar="rho",
             comp.names=c("S","I","R"),
             ic.names=c("S.0","I.0","R.0"),
             parameter.transform=function (params, logvar, logitvar, ic.names, ...) {
               params[logvar] <- exp(params[logvar])
               params[logitvar] <- plogis(params[logitvar])
               params[ic.names] <- params[ic.names]/sum(params[ic.names])
               params
             },
             parameter.inv.transform=function (params, logvar, logitvar, ic.names, ...) {
               params[ic.names] <- params[ic.names]/sum(params[ic.names])
               params[logvar] <- log(params[logvar])
               params[logitvar] <- qlogis(params[logitvar])
               params
             },
             initializer=function(params, t0, comp.names, ic.names, ...) {
               x0 <- numeric(5)
               names(x0) <- c("S","I","R","cases","W")
               fracs <- params[ic.names]
               x0[comp.names] <- round(params['pop']*fracs/sum(fracs))
               x0
             }
             )

  coef(po) <- c(                                                      
                gamma=26,mu=0.02,iota=0.01,                           
                beta1=1200,beta2=1800,beta3=600,                      
                beta.sd=1e-3,                                         
                pop=2.1e6,                                            
                rho=0.6,                                              
                S.0=26/1200,I.0=0.001,R.0=1-26/1200,             
                nbasis=3,degree=3,period=1
                )

  dyn.load(solib) ## load the shared-object library

  ## compute a trajectory of the deterministic skeleton
  tic <- Sys.time()
  X <- trajectory(po,hmax=1/52,as.data.frame=TRUE)
  toc <- Sys.time()
  print(toc-tic)

  plot(cases~time,data=X,type='l')

  ## simulate from the model
  tic <- Sys.time()
  x <- simulate(po,nsim=3,as.data.frame=TRUE)
  toc <- Sys.time()
  print(toc-tic)
  
  plot(cases~time,data=x,col=as.factor(x$sim),pch=16)

  dyn.unload(solib)

}
