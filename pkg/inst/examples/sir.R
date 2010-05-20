require(pomp)

## Coding up the SIR example using native routines results in much faster computations.
## The C codes "sir_euler_simulator", "sir_ODE", "binom_rmeasure", and "binom_dmeasure"
## are included in the "examples" directory (file "sir.c")

if (Sys.info()['sysname']=='Linux') {   # only run this under linux

  params <- c(                                                      
              gamma=26,mu=0.02,iota=0.01,                           
              beta1=1200,beta2=1800,beta3=600,                      
              beta.sd=1e-3,                                         
              pop=2.1e6,                                            
              rho=0.6,                                              
              S.0=26/1200,I.0=0.001,R.0=1-0.001-26/1200             
              )                                                     


  model <- "sir"
  pkg <- "pomp"
  modelfile <- paste(model,".c",sep="")
  cppflags <- paste("PKG_CPPFLAGS=-I",system.file("include",package=pkg),sep="")
  pkglib <- system.file(paste("libs/",pkg,.Platform$dynlib.ext,sep=""),package=pkg)
  solib <- paste(model,.Platform$dynlib.ext,sep="")

  ## compile the model into a shared-object library
  if (!file.copy(from=system.file(paste("examples/",modelfile,sep=""),package=pkg),to=getwd()))
    stop("cannot copy source code ",modelfile," to ",getwd())
  cmd <- paste(cppflags,R.home("bin/R"),"CMD SHLIB -o",solib,modelfile,pkglib)
  rv <- system(cmd)
  if (rv!=0)
    stop("cannot compile shared-object library ",solib)

  po <- pomp(
             data=data.frame(
               time=seq(from=1/52,to=4,by=1/52),
               reports=NA
               ),
             times="time",
             t0=0,
             rprocess=euler.sim(
               step.fun="sir_euler_simulator", # native routine for the simulation step
               delta.t=1/52/20
               ),
             skeleton.vectorfield="sir_ODE", # native routine for the skeleton
             rmeasure="sir_binom_rmeasure", # binomial measurement model
             dmeasure="sir_binom_dmeasure", # binomial measurement model
             PACKAGE="sir", ## name of the shared-object library
             ## the order of the observables assumed in the native routines:
             obsnames="reports",
             ## the order of the state variables assumed in the native routines:
             statenames=c("S","I","R","cases","W"),
             ## the order of the parameters assumed in the native routines:
             paramnames=c("gamma","mu","iota","beta1","beta.sd","pop","rho","nbasis","degree","period"), 
             ## reset cases to zero after each new observation:
             zeronames=c("cases"),      
             initializer=function(params,t0,...){
               p <- exp(params)
               with(
                    as.list(p),
                    {
                      fracs <- c(S.0,I.0,R.0)
                      x0 <- c(
                              round(pop*fracs/sum(fracs)), # make sure the three compartments sum to 'pop' initially
                              rep(0,2)	# zeros for 'cases' and 'W'
                              )
                      names(x0) <- c("S","I","R","cases","W")
                      x0
                    }
                    )
             }
             )

  dyn.load(solib) ## load the shared-object library

  ## compute a trajectory of the deterministic skeleton
  tic <- Sys.time()
  X <- trajectory(po,c(log(params),nbasis=3,degree=3,period=1),hmax=1/52)
  toc <- Sys.time()
  print(toc-tic)

  ## simulate from the model
  tic <- Sys.time()
  x <- simulate(po,params=c(log(params),nbasis=3,degree=3,period=1),nsim=3)
  toc <- Sys.time()
  print(toc-tic)
  
  dyn.unload(solib)

}
