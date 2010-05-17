require(pomp)

## basis functions for the seasonality
tbasis <- seq(0,25,by=1/52)
basis <- periodic.bspline.basis(tbasis,nbasis=3)
colnames(basis) <- paste("seas",1:3,sep='')

## some parameters
params <- c(
            gamma=26,mu=0.02,iota=0.01,
            beta1=1200,beta2=1800,beta3=600,
            beta.sd=1e-3,
            pop=2.1e6,
            rho=0.6,
            S.0=26/1200,I.0=0.001,R.0=1-0.001-26/1200
            )

## set up the pomp object
po <- pomp(
           times=1/52*seq.int(length=4*52),
           data=rbind(measles=numeric(52*4)),
           t0=0,
           tcovar=tbasis,
           covar=basis,
           zeronames=c("cases"),
           rprocess=euler.sim(
             step.fun=function(t,x,params,covars,delta.t,...) {
               params <- exp(params)
               with(
                    as.list(c(x,params)),
                    {
                      beta <- exp(sum(log(c(beta1,beta2,beta3))*covars))
                      beta.var <- beta.sd^2
                      dW <- rgamma(n=1,shape=delta.t/beta.var,scale=beta.var)
                      foi <- (iota+beta*I*dW/delta.t)/pop
                      trans <- c(
                                 rpois(n=1,lambda=mu*pop*delta.t),
                                 reulermultinom(n=1,size=S,rate=c(foi,mu),dt=delta.t),
                                 reulermultinom(n=1,size=I,rate=c(gamma,mu),dt=delta.t),
                                 reulermultinom(n=1,size=R,rate=c(mu),dt=delta.t)
                                 )
                      c(
                        S=S+trans[1]-trans[2]-trans[3],
                        I=I+trans[2]-trans[4]-trans[5],
                        R=R+trans[4]-trans[6],
                        cases=cases+trans[4],
                        W=if (beta.sd>0) W+(dW-delta.t)/beta.sd else W
                        )
                    }
                    )
             },
             delta.t=1/52/20
             ),
           skeleton.vectorfield=function(x,t,params,covars,...) {
             params <- exp(params)
             with(
                  as.list(c(x,params)),
                  {
                    beta <- exp(sum(log(c(beta1,beta2,beta3))*covars))
                    foi <- (iota+beta*I)/pop
                    terms <- c(
                               mu*pop,
                               foi*S,
                               mu*S,
                               gamma*I,
                               mu*I,
                               mu*R
                               )
                    xdot <- c(
                              terms[1]-terms[2]-terms[3],
                              terms[2]-terms[4]-terms[5],
                              terms[4]-terms[6],
                              terms[4]
                              )
                    ifelse(x>0,xdot,0)
                  }
                  )
           },
           measurement.model=measles~binom(size=cases,prob=exp(rho)),
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

## alternatively, one can define the computationally intensive bits using native routines:
## the C codes "sir_euler_simulator", "sir_ODE", "binom_rmeasure", and "binom_dmeasure"
## are included in the "examples" directory (file "sir.c")

if (Sys.info()['sysname']=='Linux') {

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
             times=seq(1/52,4,by=1/52),
             data=rbind(measles=numeric(52*4)),
             t0=0,
             statenames=c("S","I","R","cases","W"),
             paramnames=c("gamma","mu","iota","beta1","beta.sd","pop","rho","nbasis","degree","period"),
             zeronames=c("cases"),
             rprocess=euler.sim(
               step.fun="sir_euler_simulator",
               delta.t=1/52/20
               ),
             skeleton.vectorfield="sir_ODE",
             rmeasure="binom_rmeasure",
             dmeasure="binom_dmeasure",
             PACKAGE="sir", ## name of the shared-object library
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
