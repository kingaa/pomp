### R code from vignette source 'advanced_topics_in_pomp.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: set-opts
###################################################
  glop <- options(keep.source=TRUE,width=60,continue=" ",prompt=" ")
  library(pomp)
  pdf.options(useDingbats=FALSE)
  set.seed(5384959)


###################################################
### code chunk number 2: plugin-R-code
###################################################
pompExample(ou2)
ou2.dat <- as.data.frame(ou2)

pomp( 
     data=ou2.dat[c("time","y1","y2")],
     times="time",
     t0=0,
     rprocess=discrete.time.sim(
       step.fun=function (x, t, params, ...) {
         eps <- rnorm(n=2,mean=0,sd=1) # noise terms
         xnew <- c(
                   x1=params["alpha.1"]*x["x1"]+params["alpha.3"]*x["x2"]+
                       params["sigma.1"]*eps[1],
                   x2=params["alpha.2"]*x["x1"]+params["alpha.4"]*x["x2"]+
                       params["sigma.2"]*eps[1]+params["sigma.3"]*eps[2]
                   )
         names(xnew) <- c("x1","x2")
         xnew
       }
       )
     ) -> ou2.Rplug


###################################################
### code chunk number 3: plugin-R-code-sim (eval = FALSE)
###################################################
## simdat.Rplug <- simulate(ou2.Rplug,params=coef(ou2),nsim=5000,states=T)


###################################################
### code chunk number 4: plugin-R-code-eval
###################################################
binary.file <- "plugin-R-code.rda"
if (file.exists(binary.file)) {
  load(binary.file)
} else {
  tic <- Sys.time()
simdat.Rplug <- simulate(ou2.Rplug,params=coef(ou2),nsim=5000,states=T)
  toc <- Sys.time()
  etime.Rplug <- toc-tic
  n.Rplug <- dim(simdat.Rplug)[2]
  save(etime.Rplug,n.Rplug,file=binary.file,compress='xz')
}


###################################################
### code chunk number 5: vectorized-R-code (eval = FALSE)
###################################################
## ou2.Rvect.rprocess <- function (xstart, times, params, ...) {
##   nrep <- ncol(xstart)                  # number of realizations
##   ntimes <- length(times)               # number of timepoints
##   ## unpack the parameters (for legibility only)
##   alpha.1 <- params["alpha.1",]
##   alpha.2 <- params["alpha.2",]
##   alpha.3 <- params["alpha.3",]
##   alpha.4 <- params["alpha.4",]
##   sigma.1 <- params["sigma.1",]
##   sigma.2 <- params["sigma.2",]
##   sigma.3 <- params["sigma.3",]
##   ## x is the array of states to be returned: it must have rownames
##   x <- array(0,dim=c(2,nrep,ntimes))
##   rownames(x) <- rownames(xstart)
##   ## xnow holds the current state values
##   x[,,1] <- xnow <- xstart
##   tnow <- times[1]
##   for (k in seq.int(from=2,to=ntimes,by=1)) {
##     tgoal <- times[k]
##     while (tnow < tgoal) {              # take one step at a time
##       eps <- array(rnorm(n=2*nrep,mean=0,sd=1),dim=c(2,nrep))
##       tmp <- alpha.1*xnow['x1',]+alpha.3*xnow['x2',]+
##         sigma.1*eps[1,]
##       xnow['x2',] <- alpha.2*xnow['x1',]+alpha.4*xnow['x2',]+
##         sigma.2*eps[1,]+sigma.3*eps[2,]
##       xnow['x1',] <- tmp
##       tnow <- tnow+1
##     }
##     x[,,k] <- xnow
##   }
##   x
## }


###################################################
### code chunk number 6: vectorized-R-pomp (eval = FALSE)
###################################################
## ou2.Rvect <- pomp(ou2.Rplug,rprocess=ou2.Rvect.rprocess)


###################################################
### code chunk number 7: theta (eval = FALSE)
###################################################
## theta <- c(
##            x1.0=-3, x2.0=4,
##            tau=1,
##            alpha.1=0.8, alpha.2=-0.5, alpha.3=0.3, alpha.4=0.9,
##            sigma.1=3, sigma.2=-0.5, sigma.3=2
##            )


###################################################
### code chunk number 8: vectorized-R-code-sim (eval = FALSE)
###################################################
## simdat.Rvect <- simulate(ou2.Rvect,params=theta,states=T,nsim=100000)


###################################################
### code chunk number 9: vectorized-R-code-eval
###################################################
binary.file <- "vectorized-R-code.rda"
if (file.exists(binary.file)) {
  load(binary.file)
} else {
ou2.Rvect.rprocess <- function (xstart, times, params, ...) {
  nrep <- ncol(xstart)                  # number of realizations
  ntimes <- length(times)               # number of timepoints
  ## unpack the parameters (for legibility only)
  alpha.1 <- params["alpha.1",]
  alpha.2 <- params["alpha.2",]
  alpha.3 <- params["alpha.3",]
  alpha.4 <- params["alpha.4",]
  sigma.1 <- params["sigma.1",]
  sigma.2 <- params["sigma.2",]
  sigma.3 <- params["sigma.3",]
  ## x is the array of states to be returned: it must have rownames
  x <- array(0,dim=c(2,nrep,ntimes))
  rownames(x) <- rownames(xstart)
  ## xnow holds the current state values
  x[,,1] <- xnow <- xstart
  tnow <- times[1]
  for (k in seq.int(from=2,to=ntimes,by=1)) {
    tgoal <- times[k]
    while (tnow < tgoal) {              # take one step at a time
      eps <- array(rnorm(n=2*nrep,mean=0,sd=1),dim=c(2,nrep))
      tmp <- alpha.1*xnow['x1',]+alpha.3*xnow['x2',]+
        sigma.1*eps[1,]
      xnow['x2',] <- alpha.2*xnow['x1',]+alpha.4*xnow['x2',]+
        sigma.2*eps[1,]+sigma.3*eps[2,]
      xnow['x1',] <- tmp
      tnow <- tnow+1
    }
    x[,,k] <- xnow
  }
  x
}
ou2.Rvect <- pomp(ou2.Rplug,rprocess=ou2.Rvect.rprocess)
theta <- c(
           x1.0=-3, x2.0=4,
           tau=1,
           alpha.1=0.8, alpha.2=-0.5, alpha.3=0.3, alpha.4=0.9,
           sigma.1=3, sigma.2=-0.5, sigma.3=2
           )
  tic <- Sys.time()
simdat.Rvect <- simulate(ou2.Rvect,params=theta,states=T,nsim=100000)
  toc <- Sys.time()
etime.Rvect <- toc-tic
units(etime.Rvect) <- units(etime.Rplug)
n.Rvect <- dim(simdat.Rvect)[2]
save(etime.Rvect,n.Rvect,file=binary.file,compress='xz')
}


###################################################
### code chunk number 10: view-ou2-source (eval = FALSE)
###################################################
## file.show(file=system.file("examples/ou2.c",package="pomp"))


###################################################
### code chunk number 11: view-pomp.h (eval = FALSE)
###################################################
## file.show(file=system.file("include/pomp.h",package="pomp"))


###################################################
### code chunk number 12: plugin-C-code (eval = FALSE)
###################################################
## ou2.Cplug <- pomp(
##                   ou2.Rplug,
##                   rprocess=discrete.time.sim("ou2_step"),
##                   paramnames=c(
##                     "alpha.1","alpha.2","alpha.3","alpha.4",
##                     "sigma.1","sigma.2","sigma.3",
##                     "tau"
##                     ),
##                   statenames=c("x1","x2"),
##                   obsnames=c("y1","y2")
##                   )


###################################################
### code chunk number 13: plugin-C-sim (eval = FALSE)
###################################################
## simdat.Cplug <- simulate(ou2.Cplug,params=theta,states=T,nsim=100000)


###################################################
### code chunk number 14: plugin-C-sim-eval
###################################################
binary.file <- "plugin-C-code.rda"
if (file.exists(binary.file)) {
  load(binary.file)
} else {
ou2.Cplug <- pomp(
                  ou2.Rplug,
                  rprocess=discrete.time.sim("ou2_step"),
                  paramnames=c(
                    "alpha.1","alpha.2","alpha.3","alpha.4",
                    "sigma.1","sigma.2","sigma.3",
                    "tau"
                    ),
                  statenames=c("x1","x2"),
                  obsnames=c("y1","y2")
                  )
  tic <- Sys.time()
simdat.Cplug <- simulate(ou2.Cplug,params=theta,states=T,nsim=100000)
  toc <- Sys.time()
etime.Cplug <- toc-tic
  n.Cplug <- dim(simdat.Cplug)[2]
units(etime.Cplug) <- units(etime.Rplug)
speedup <- as.numeric(etime.Rplug/n.Rplug)/as.numeric(etime.Cplug/n.Cplug)
save(n.Cplug,etime.Cplug,speedup,file=binary.file,compress='xz')
}


###################################################
### code chunk number 15: advanced_topics_in_pomp.Rnw:282-283 (eval = FALSE)
###################################################
## file.show(file=system.file("examples/ou2.c",package="pomp"))


###################################################
### code chunk number 16: vectorized-C-code
###################################################
ou2.Cvect.rprocess <- function (xstart, times, params, ...) {
  nvar <- nrow(xstart)
  npar <- nrow(params)
  nrep <- ncol(xstart)
  ntimes <- length(times)
  array(
        .C("ou2_adv",
           X = double(nvar*nrep*ntimes),
           xstart = as.double(xstart),
           par = as.double(params),
           times = as.double(times),
           n = as.integer(c(nvar,npar,nrep,ntimes))
           )$X,
        dim=c(nvar,nrep,ntimes),
        dimnames=list(rownames(xstart),NULL,NULL)        
        )
}


###################################################
### code chunk number 17: vectorized-C-code-pomp
###################################################
ou2.Cvect <- pomp(
                  ou2.Rplug,
                  rprocess=ou2.Cvect.rprocess
                  )
paramnames <- c(
                "alpha.1","alpha.2","alpha.3","alpha.4",
                "sigma.1","sigma.2","sigma.3",
                "tau",
                "x1.0","x2.0"
                )


###################################################
### code chunk number 18: vectorized-C-code-sim (eval = FALSE)
###################################################
## simdat.Cvect <- simulate(ou2.Cvect,params=theta[paramnames],nsim=100000,states=T)


###################################################
### code chunk number 19: vectorized-C-code-eval
###################################################
binary.file <- "vectorized-C-code.rda"
if (file.exists(binary.file)) {
  load(binary.file)
} else {
  tic <- Sys.time()
simdat.Cvect <- simulate(ou2.Cvect,params=theta[paramnames],nsim=100000,states=T)
  toc <- Sys.time()
etime.Cvect <- toc-tic
  n.Cvect <- dim(simdat.Cvect)[2]
units(etime.Cvect) <- units(etime.Rplug)
speedup <- as.numeric(etime.Rplug/n.Rplug)/as.numeric(etime.Cvect/n.Cvect)
save(n.Cvect,etime.Cvect,speedup,file=binary.file,compress='xz')
}


###################################################
### code chunk number 20: pomp-builder-measmod
###################################################

## negative binomial measurement model
## E[cases|incid] = rho*incid
## Var[cases|incid] = rho*incid*(1+rho*incid/theta)
rmeas <- '
  cases = rnbinom_mu(theta,rho*incid);
'

dmeas <- '
  lik = dnbinom_mu(cases,theta,rho*incid,give_log);
'



###################################################
### code chunk number 21: pomp-builder-stepfn
###################################################

## SIR process model with extra-demographic stochasticity
## and seasonal transmission
step.fn <- '
  int nrate = 6;
  double rate[nrate];		// transition rates
  double trans[nrate];		// transition numbers
  double beta;			// transmission rate
  double dW;			// white noise increment
  int k;

  // seasonality in transmission
  beta = beta1*seas1+beta2*seas2+beta3*seas3;

  // compute the environmental stochasticity
  dW = rgammawn(beta_sd,dt);

  // compute the transition rates
  rate[0] = mu*popsize;		// birth into susceptible class
  rate[1] = (iota+beta*I*dW/dt)/popsize; // force of infection
  rate[2] = mu;			// death from susceptible class
  rate[3] = gamma;		// recovery
  rate[4] = mu;			// death from infectious class
  rate[5] = mu; 		// death from recovered class

  // compute the transition numbers
  trans[0] = rpois(rate[0]*dt);	// births are Poisson
  reulermultinom(2,S,&rate[1],dt,&trans[1]);
  reulermultinom(2,I,&rate[3],dt,&trans[3]);
  reulermultinom(1,R,&rate[5],dt,&trans[5]);

  // balance the equations
  S += trans[0]-trans[1]-trans[2];
  I += trans[1]-trans[3]-trans[4];
  R += trans[3]-trans[5];
  incid += trans[3];		// incidence is cumulative recoveries
  if (beta_sd > 0.0) W += (dW-dt)/beta_sd; // increment has mean = 0, variance = dt
'



###################################################
### code chunk number 22: pomp-builder-skel
###################################################

skel <- '
  int nrate = 6;
  double rate[nrate];		// transition rates
  double term[nrate];		// transition numbers
  double beta;			// transmission rate
  double dW;			// white noise increment
  int k;
  
  beta = beta1*seas1+beta2*seas2+beta3*seas3;

  // compute the transition rates
  rate[0] = mu*popsize;		// birth into susceptible class
  rate[1] = (iota+beta*I)/popsize; // force of infection
  rate[2] = mu;			// death from susceptible class
  rate[3] = gamma;		// recovery
  rate[4] = mu;			// death from infectious class
  rate[5] = mu; 		// death from recovered class

  // compute the several terms
  term[0] = rate[0];
  term[1] = rate[1]*S;
  term[2] = rate[2]*S;
  term[3] = rate[3]*I;
  term[4] = rate[4]*I;
  term[5] = rate[5]*R;

  // assemble the differential equations
  DS = term[0]-term[1]-term[2];
  DI = term[1]-term[3]-term[4];
  DR = term[3]-term[5];
  Dincid = term[3];		// accumulate the new I->R transitions
  DW = 0;
'



###################################################
### code chunk number 23: pomp-builder-partrans
###################################################

## parameter transformations
## note we use barycentric coordinates for the initial conditions
## the success of this depends on S0, I0, R0 being in
## adjacent memory locations, in that order
partrans <- "
  Tgamma = exp(gamma);
  Tmu = exp(mu);
  Tiota = exp(iota);
  Tbeta1 = exp(beta1);
  Tbeta2 = exp(beta2);
  Tbeta3 = exp(beta3);
  Tbeta_sd = exp(beta_sd);
  Trho = expit(rho);
  Ttheta = exp(theta);
  from_log_barycentric(&TS_0,&S_0,3);
"

paruntrans <- "
  Tgamma = log(gamma);
  Tmu = log(mu);
  Tiota = log(iota);
  Tbeta1 = log(beta1);
  Tbeta2 = log(beta2);
  Tbeta3 = log(beta3);
  Tbeta_sd = log(beta_sd);
  Trho = logit(rho);
  Ttheta = log(theta);
  to_log_barycentric(&TS_0,&S_0,3);
"



###################################################
### code chunk number 24: pomp-builder-covar
###################################################

covartab <- data.frame(
                       time=seq(from=-1/52,to=10+1/52,by=1/26)
                       )

covartab <- cbind(
                  covartab,
                  with(covartab,
                       periodic.bspline.basis(
                                              x=time,
                                              nbasis=3,
                                              degree=3,
                                              period=1,
                                              names="seas%d"
                                              )
                       )
                  )



###################################################
### code chunk number 25: pomp-builder (eval = FALSE)
###################################################
## 
## pompBuilder(
##             name="SIR",
##             data=data.frame(
##               cases=NA,
##               time=seq(0,10,by=1/52)
##               ),
##             times="time",
##             t0=-1/52,
##             dmeasure=dmeas,
##             rmeasure=rmeas,
##             step.fn=step.fn,
##             step.fn.delta.t=1/52/20,
##             skeleton.type="vectorfield",
##             skeleton=skel,
##             covar=covartab,
##             tcovar="time",
##             parameter.transform=partrans,
##             parameter.inv.transform=paruntrans,
##             statenames=c("S","I","R","incid","W"),
##             paramnames=c(
##               "gamma","mu","iota",
##               "beta1","beta2","beta3","beta.sd",
##               "popsize","rho","theta",
##               "S.0","I.0","R.0"
##               ), 
##             zeronames=c("incid","W"),
##             initializer=function(params, t0, ...) {
##               x0 <- setNames(numeric(5),c("S","I","R","incid","W"))
##               fracs <- params[c("S.0","I.0","R.0")]
##               x0[1:3] <- round(params['popsize']*fracs/sum(fracs))
##               x0
##             }
##             ) -> sir
## 


###################################################
### code chunk number 26: pomp-builder-eval
###################################################
if (Sys.getenv("POMP_BUILD_VIGNETTES")=="yes") {
  require(pomp)

pompBuilder(
            name="SIR",
            data=data.frame(
              cases=NA,
              time=seq(0,10,by=1/52)
              ),
            times="time",
            t0=-1/52,
            dmeasure=dmeas,
            rmeasure=rmeas,
            step.fn=step.fn,
            step.fn.delta.t=1/52/20,
            skeleton.type="vectorfield",
            skeleton=skel,
            covar=covartab,
            tcovar="time",
            parameter.transform=partrans,
            parameter.inv.transform=paruntrans,
            statenames=c("S","I","R","incid","W"),
            paramnames=c(
              "gamma","mu","iota",
              "beta1","beta2","beta3","beta.sd",
              "popsize","rho","theta",
              "S.0","I.0","R.0"
              ), 
            zeronames=c("incid","W"),
            initializer=function(params, t0, ...) {
              x0 <- setNames(numeric(5),c("S","I","R","incid","W"))
              fracs <- params[c("S.0","I.0","R.0")]
              x0[1:3] <- round(params['popsize']*fracs/sum(fracs))
              x0
            }
            ) -> sir

}


###################################################
### code chunk number 27: sir-sim (eval = FALSE)
###################################################
## 
## coef(sir) <- c(
##                gamma=26,mu=0.02,iota=0.01,
##                beta1=400,beta2=480,beta3=320,
##                beta.sd=0.001,
##                popsize=2.1e6,
##                rho=0.6,theta=10,
##                S.0=26/400,I.0=0.001,R.0=1
##                )
## 
## sir <- simulate(sir,seed=3493885L)
## traj <- trajectory(sir,hmax=1/52)
## 


###################################################
### code chunk number 28: sir-sim-eval
###################################################
binary.file <- "sim-sim.rda"
if (file.exists(binary.file)) {
  load(binary.file)
} else {

coef(sir) <- c(
               gamma=26,mu=0.02,iota=0.01,
               beta1=400,beta2=480,beta3=320,
               beta.sd=0.001,
               popsize=2.1e6,
               rho=0.6,theta=10,
               S.0=26/400,I.0=0.001,R.0=1
               )

sir <- simulate(sir,seed=3493885L)
traj <- trajectory(sir,hmax=1/52)

  save(sir,traj,file=binary.file,compress='xz')
}


###################################################
### code chunk number 29: sir-plot
###################################################
plot(sir)


###################################################
### code chunk number 30: advanced_topics_in_pomp.Rnw:631-638
###################################################
pompExample(ou2)
true.p <- coef(ou2)
x0 <- init.state(ou2)
x0
new.p <- cbind(true.p,true.p,true.p)
new.p["x1.0",] <- 1:3
init.state(ou2,params=new.p)


###################################################
### code chunk number 31: advanced_topics_in_pomp.Rnw:647-650
###################################################
x <- rprocess(ou2,xstart=x0,times=time(ou2,t0=T),params=true.p)
dim(x)
x[,,1:5]


###################################################
### code chunk number 32: advanced_topics_in_pomp.Rnw:658-662
###################################################
x <- x[,,-1,drop=F]
y <- rmeasure(ou2,x=x,times=time(ou2),params=true.p)
dim(y)
y[,,1:5]


###################################################
### code chunk number 33: advanced_topics_in_pomp.Rnw:668-671
###################################################
fp <- dprocess(ou2,x=x,times=time(ou2),params=true.p)
dim(fp)
fp[,36:40]


###################################################
### code chunk number 34: advanced_topics_in_pomp.Rnw:673-676
###################################################
fm <- dmeasure(ou2,y=y[,1,],x=x,times=time(ou2),params=true.p)
dim(fm)
fm[,36:40]


###################################################
### code chunk number 35: all-examples (eval = FALSE)
###################################################
## pompExample()


###################################################
### code chunk number 36: restore-opts
###################################################
options(glop)


