library(pomp)

modelfile <- system.file("examples/sir.c",package="pomp")
includedir <- system.file("include",package="pomp")
lib <- system.file("libs/pomp.so",package="pomp")

## compile the model into shared-object library
system(paste("gcc -fPIC -c -I/usr/local/lib64/R/include -I",includedir,modelfile))
system(paste("R CMD SHLIB -o sir.so sir.o",lib))

## set up a lookup table for basis functions for the seasonality
tbasis <- seq(0,20,by=1/52)
basis <- periodic.bspline.basis(tbasis,nbasis=3)

## some parameters
params <- c(gamma=26,mu=0.2,iota=0.01,
            beta1=1200,beta2=2100,beta3=300,
            beta.sd=0.1,
            pop=2.1e5,
            rho=0.6,
            S.0=26/1200,I.0=0.001,R.0=1-0.001-26/1200
            )

dyn.load("sir.so")                    # load the shared-object library

## set up the pomp object
po <- pomp(
           time=seq(1/52,20,by=1/52),
           data=rbind(measles=numeric(52*20)),
           t0=0,
           tbasis=tbasis,
           basis=basis,
           dt=1/52/20,
           rprocess=function(xstart,times,params,dt,tbasis,basis,...){
             euler.simulate(
                            xstart=xstart,
                            times=times,
                            params=params,
                            euler.step.fun="sir_euler_multinomial",
                            delta.t=dt,
                            statenames=c("S","I","R","cases","W","trans1"),
                            paramnames=c("gamma","mu","iota","beta1","beta.sd","pop"),
                            zeronames=c("cases"),
                            tcovar=tbasis,
                            covar=basis,
                            PACKAGE="sir"
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
             array(
                   data=dbinom(
                     x=y,
                     size=x["cases",,],
                     prob=exp(params["rho",]),
                     log=log
                     ),
                   dim=c(nreps,ntimes)
                   )
           },
           initializer=function(params,t0,...){
             p <- exp(params)
             pop <- p["pop"]
             fracs <- p[c("S.0","I.0","R.0")]
             x0 <- c(
                     round(pop*fracs/sum(fracs)),
                     0,rep(0,7)
                     )
             names(x0) <- c(
                            "S","I","R","cases","W",
                            paste("trans",1:6,sep="")
                            )
             x0
           }
           )

## simulate from the model
tic <- Sys.time()
x <- simulate(po,params=log(params),nsim=10)
toc <- Sys.time()
print(toc-tic)

X <- simulate(po,params=log(params),nsim=10,states=T)
nm <- rownames(X)
dim(X) <- c(11,10,1041)
rownames(X) <- nm

f <- dmeasure(
              po,
              y=data.array(x[[1]]),
              x=X[,,-1,drop=F],
              times=time(po),
              params=matrix(
                log(params),
                nrow=length(params),
                ncol=10,
                dimnames=list(names(params),NULL)
                ),
              log=T
              )


dyn.unload("sir.so")
