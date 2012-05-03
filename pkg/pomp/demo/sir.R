require(pomp)

dmeas <- "
  lik = dbinom(reports,nearbyint(cases),rho,give_log);
"
rmeas <- "
  reports = rbinom(nearbyint(cases),rho);
"
step.fn <- '
  int nrate = 6;
  double rate[nrate];		// transition rates
  double trans[nrate];		// transition numbers
  double beta;			// transmission rate
  double dW;			// white noise increment
  double period = 1.0;          // period of the seasonality
  int nbasis = 3;               // number of seasonality basis functions
  int deg = 3;                  // degree of the B-spline basis functions
  double seasonality[nbasis];
  int k;

  // compute transmission rate from seasonality
  periodic_bspline_basis_eval(t,period,deg,nbasis,&seasonality[0]); // evaluate the periodic B-spline basis
  for (k = 0, beta = 0; k < nbasis; k++) 
    beta += (&beta1)[k]*seasonality[k];

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
  cases += trans[3];		// cases are cumulative recoveries
  if (beta_sd > 0.0) W += (dW-dt)/beta_sd; // increment has mean = 0, variance = dt
'
skel <- '
  int nrate = 6;
  double rate[nrate];		// transition rates
  double term[nrate];		// transition numbers
  double beta;			// transmission rate
  double dW;			// white noise increment
  double period = 1.0;          // period of the seasonality
  int nbasis = 3;
  int deg = 3;
  double seasonality[nbasis];
  int k;
  
  // compute transmission rate from seasonality
  if (nbasis <= 0) return;
  periodic_bspline_basis_eval(t,period,deg,nbasis,&seasonality[0]); // evaluate the periodic B-spline basis
  for (k = 0, beta = 0; k < nbasis; k++) 
    beta += (&beta1)[k]*seasonality[k];

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
  Dcases = term[3];		// accumulate the new I->R transitions
  DW = 0;
'

pompBuilder(
            name="SIR",
            data=data.frame(
              time=seq(from=1/52,to=4,by=1/52),
              reports=NA
              ),
            times="time",
            t0=0,
            dmeasure=dmeas,
            rmeasure=rmeas,
            step.fn=step.fn,
            step.fn.delta.t=1/52/20,
            skeleton.type="vectorfield",
            skeleton=skel,
            statenames=c("S","I","R","cases","W"),
            paramnames=c(
              "gamma","mu","iota","beta1","beta.sd","popsize","rho"
              ), 
            zeronames=c("cases","W"),
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
              x0[comp.names] <- round(params['popsize']*fracs/sum(fracs))
              x0
            }
            ) -> po

coef(po) <- c(
              gamma=26,mu=0.02,iota=0.01,
              beta1=120,beta2=140,beta3=100,
              beta.sd=1e-3,
              popsize=5e5,
              rho=0.6,
              S.0=26/120,I.0=0.001,R.0=1-26/120
              )

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

dyn.unload("SIR.so")
