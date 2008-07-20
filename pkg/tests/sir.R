library(pomp)

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
## the C codes "sir_euler_simulator" and "sir_euler_density" are included in the "examples" directory (file "sir.c")
po <- pomp(
           times=1/52*seq.int(length=4*52),
           data=rbind(measles=numeric(52*4)),
           t0=0,
           tcovar=tbasis,
           covar=basis,
           delta.t=1/52/20,
           zeronames=c("cases"),
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
                      W=if (beta.sd>0) W+(dW-delta.t)/beta.sd else W,
                      B=trans[1],
                      SI=trans[2],
                      SD=trans[3],
                      IR=trans[4],
                      ID=trans[5],
                      RD=trans[6],
                      dW=dW
                      )
                  }
                  )
           },
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
           rprocess=euler.simulate,
           dprocess=euler.density,
           measurement.model=measles~binom(size=cases,prob=exp(rho)),
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

t <- seq(0,4/52,by=1/52/25)
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
print(apply(f,1,sum),digits=4)

g <- dmeasure(
              po,
              y=rbind(measles=X$obs[,7,]),
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
print(apply(g,1,sum),digits=4)

t <- seq(0,2,by=1/52)
X <- simulate(po,params=log(params),nsim=1,states=TRUE,obs=TRUE,times=t)

h <- skeleton(
              po,
              x=X$states[,1,55:70,drop=FALSE],
              t=t[55:70],
              params=as.matrix(log(params))
              )
print(h[c("S","I","R"),,],digits=4)

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
           step.fun="sir_euler_simulator",
           rprocess=euler.simulate,
           dens.fun="sir_euler_density",
           dprocess=euler.density,
           skeleton="sir_ODE",
           PACKAGE="pomp",
           measurement.model=measles~binom(size=cases,prob=exp(rho)),
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

t <- seq(0,4/52,by=1/52/25)
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
print(apply(f,1,sum),digits=4)

g <- dmeasure(
              po,
              y=rbind(measles=X$obs[,7,]),
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
print(apply(g,1,sum),digits=4)

t <- seq(0,2,by=1/52)
X <- simulate(po,params=log(params),nsim=1,states=TRUE,obs=TRUE,times=t)

h <- skeleton(
              po,
              x=X$states[,1,55:70,drop=FALSE],
              t=t[55:70],
              params=as.matrix(log(params))
              )
print(h[c("S","I","R"),,],digits=4)
