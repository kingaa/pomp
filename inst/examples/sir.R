require(pomp)

modelfile <- system.file("examples/sir.c",package="pomp")
includedir <- system.file("include",package="pomp")
lib <- system.file("libs/pomp.so",package="pomp")

## compile the model into shared-object library
system(paste("cp",modelfile,"."))
system(paste("cp ",includedir,"/euler.h .",sep=""))
system(paste("cp ",includedir,"/lookup_table.h .",sep=""))
system(paste("R CMD SHLIB -o sir.so sir.c",lib))

## set up a lookup table for basis functions for the seasonality
tbasis <- seq(0,20,by=1/52)
basis <- periodic.bspline.basis(tbasis,nbasis=3)
colnames(basis) <- paste("seas",1:3,sep='.')

## some parameters
params <- c(
            gamma=26,mu=0.2,iota=0.01,
            beta1=1200,beta2=2100,beta3=300,
            beta.sd=0.1,
            pop=2.1e5,
            rho=0.6,
            S.0=26/1200,I.0=0.001,R.0=1-0.001-26/1200
            )

## set up the pomp object
po <- pomp(
           time=seq(1/52,20,by=1/52),
           data=rbind(measles=numeric(52*20)),
           t0=0,
           tbasis=tbasis,
           basis=basis,
           dt=1/52/20,
           statenames=c("S","I","R","cases","W","trans1"),
           paramnames=c("gamma","mu","iota","beta1","beta.sd","pop"),
           basisnames=c("seas.1"),
           zeronames=c("cases"),
           rprocess=function(xstart,times,params,dt,tbasis,basis,
             statenames,paramnames,basisnames,zeronames,...){
             euler.simulate(
                            xstart=xstart,
                            times=times,
                            params=params,
                            euler.step.fun="sir_euler_simulator",
                            delta.t=dt,
                            statenames=statenames,
                            paramnames=paramnames,
                            covarnames=basisnames,
                            zeronames=zeronames,
                            tcovar=tbasis,
                            covar=basis
                            )
           },
           dprocess=function(x,times,params,log,tbasis,basis,
             statenames,paramnames,basisnames,...){
             euler.density(
                           x=x,
                           times=times,
                           params=params,
                           euler.dens.fun="sir_euler_density",
                           statenames=statenames,
                           paramnames=paramnames,
                           covarnames=basisnames,
                           tcovar=tbasis,
                           covar=basis,
                           log=log
                           )
           },
           rmeasure=function(x,times,params,...){
             nsims <- ncol(params)
             ntimes <- length(times)
             array(
                   data=rbinom(
                     n=nsims*ntimes,
                     size=x["cases",,],
                     prob=exp(params["rho",])
                     ),
                   dim=c(1,nsims,ntimes),
                   dimnames=list("measles",NULL,NULL)
                   )
           },
           dmeasure=function(y,x,times,params,log=FALSE,...){
             nreps <- ncol(params)
             ntimes <- length(times)
             f <- array(dim=c(nreps,ntimes))
             for (k in 1:nreps)
               f[k,] <- dbinom(
                               x=y,
                               size=x["cases",k,],
                               prob=exp(params["rho",k]),
                               log=log
                               )
             f
           },
           initializer=function(params,t0,...){
             p <- exp(params)
             pop <- p["pop"]
             fracs <- p[c("S.0","I.0","R.0")]
             x0 <- c(
                     round(pop*fracs/sum(fracs)), # make sure the three compartments sum to 'pop' initially
                     rep(0,9)	# zeros for 'cases', 'W', and the transition numbers
                     )
             names(x0) <- c(
                            "S","I","R","cases","W",
                            paste("trans",1:7,sep="")
                            )
             x0
           }
           )

dyn.load("sir.so")                    # load the shared-object library

## simulate from the model
tic <- Sys.time()
x <- simulate(po,params=log(params),nsim=3)
toc <- Sys.time()
print(toc-tic)
plot(x[[1]])

t <- seq(0,4/52,1/52/20)
X <- simulate(po,params=log(params),nsim=10,states=TRUE,obs=TRUE,times=t)

f <- dprocess(
              po,
              x=X$states[,,31:40],
              times=t[31:40],
              params=matrix(
                log(params),
                nrow=length(params),
                ncol=10,
                dimnames=list(names(params),NULL)
                ),
              log=TRUE
              )
apply(f,1,sum)

g <- dmeasure(
              po,
              y=t(as.matrix(X$obs[,7,])),
              x=X$states,
              times=t,
              params=matrix(
                log(params),
                nrow=length(params),
                ncol=10,
                dimnames=list(names(params),NULL)
                ),
              log=TRUE
              )
apply(g,1,sum)
}

dyn.unload("sir.so")
