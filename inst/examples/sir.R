require(pomp)

modelfile <- system.file("examples/sir.c",package="pomp")
includedir <- system.file("include",package="pomp")
lib <- system.file("libs/pomp.so",package="pomp")

## compile the model into shared-object library
system(paste("cp",modelfile,"."))
system(paste("cp ",includedir,"/euler.h .",sep=""))
system(paste("cp ",includedir,"/lookup_table.h .",sep=""))
system(paste("R CMD SHLIB -o sir.so sir.c",lib))

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
           times=seq(1/52,4,by=1/52),
           data=rbind(measles=numeric(52*4)),
           t0=0,
           tcovar=tbasis,
           covar=basis,
           delta.t=1/52/20,
           statenames=c("S","I","R","cases","W","B","dW"),
           paramnames=c("gamma","mu","iota","beta1","beta.sd","pop"),
           covarnames=c("seas1"),
           zeronames=c("cases"),
           measurement.model=measles~binom(size=cases,prob=exp(rho)),
           rprocess=euler.simulate,
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
                      S+trans[1]-trans[2]-trans[3],
                      I+trans[2]-trans[4]-trans[5],
                      R+trans[4]-trans[6],
                      cases+trans[4],
                      if (beta.sd>0) W+(dW-delta.t)/beta.sd else W,
                      trans,
                      dW
                      )
                  }
                  )
           },
           dprocess=euler.density,
           dens.fun=function(t1,t2,params,x1,x2,covars,...) {
             params <- exp(params)
             with(
                  as.list(params),
                  {
                    dt <- t2-t1
                    beta <- exp(sum(log(c(beta1,beta2,beta3))*covars))
                    beta.var <- beta.sd^2
                    dW <- x2['dW']
                    foi <- (iota+beta*x1["I"]*dW/dt)/pop
                    probs <- c(
                               dpois(x=x2["B"],lambda=mu*pop*dt,log=T),
                               deulermultinom(x=x2[c("SI","SD")],size=x1["S"],rate=c(foi,mu),dt=dt,log=T),
                               deulermultinom(x=x2[c("IR","ID")],size=x1["I"],rate=c(gamma,mu),dt=dt,log=T),
                               deulermultinom(x=x2["RD"],size=x1["R"],rate=c(mu),dt=dt,log=T),
                               dgamma(x=dW,shape=dt/beta.var,scale=beta.var,log=T)
                               )
                    sum(probs)
                  }
                  )
           },
           skeleton=function(x,t,params,covars,...) {
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
                    c(
                      terms[1]-terms[2]-terms[3],
                      terms[2]-terms[4]-terms[5],
                      terms[4]-terms[6],
                      terms[4]
                      )
                  }
                  )
           },
           initializer=function(params,t0,...){
             p <- exp(params)
             with(
                  as.list(p),
                  {
                    fracs <- c(S.0,I.0,R.0)
                    x0 <- c(
                            round(pop*fracs/sum(fracs)), # make sure the three compartments sum to 'pop' initially
                            rep(0,9)	# zeros for 'cases', 'W', and the transition numbers
                            )
                    names(x0) <- c("S","I","R","cases","W","B","SI","SD","IR","ID","RD","dW")
                    x0
                  }
                  )
           }
           )

set.seed(3049953)
## simulate from the model
tic <- Sys.time()
x <- simulate(po,params=log(params),nsim=3)
toc <- Sys.time()
print(toc-tic)

# alternatively, one can define the computationally intensive bits using native routines:
## the C codes "sir_euler_simulator" and "sir_euler_density" are included in the "examples" directory (file "sir.c")
po <- pomp(
           times=seq(1/52,4,by=1/52),
           data=rbind(measles=numeric(52*4)),
           t0=0,
           tcovar=tbasis,
           covar=basis,
           delta.t=1/52/20,
           statenames=c("S","I","R","cases","W","B","dW"),
           paramnames=c("gamma","mu","iota","beta1","beta.sd","pop"),
           covarnames=c("seas1"),
           zeronames=c("cases"),
           measurement.model=measles~binom(size=cases,prob=exp(rho)),
           rprocess=euler.simulate,
           step.fun="sir_euler_simulator",
           dprocess=euler.density,
           dens.fun="sir_euler_density",
           skeleton="sir_ODE",
           PACKAGE="pomp",
           initializer=function(params,t0,...){
             p <- exp(params)
             with(
                  as.list(p),
                  {
                    fracs <- c(S.0,I.0,R.0)
                    x0 <- c(
                            round(pop*fracs/sum(fracs)), # make sure the three compartments sum to 'pop' initially
                            rep(0,9)	# zeros for 'cases', 'W', and the transition numbers
                            )
                    names(x0) <- c("S","I","R","cases","W","B","SI","SD","IR","ID","RD","dW")
                    x0
                  }
                  )
           }
           )

dyn.load("sir.so")                    # load the shared-object library

## simulate from the model
tic <- Sys.time()
x <- simulate(po,params=log(params),nsim=3)
toc <- Sys.time()
print(toc-tic)
plot(x[[1]])

dyn.unload("sir.so")
