library(pomp)

## negative binomial measurement model
## E[cases|incid] = rho*incid
## Var[cases|incid] = rho*incid*(1+rho*incid/theta)
rmeas <- "
  cases = rnbinom_mu(theta,rho*incid);
"

## three basis functions
globals <- "
  static int nbasis = 3;
"

## SIR process model with extra-demographic stochasticity
## and seasonal transmission
step.fun <- "
  double rate[6];		// transition rates
  double trans[6];	// transition numbers
  double beta;			// transmission rate
  double dW;		    // white noise increment
  int k;

  // seasonality in transmission
  for (k = 0, beta = 0.0; k < nbasis; k++)
     beta += (&beta1)[k]*(&seas_1)[k];

  // compute the environmental stochasticity
  dW = rgammawn(beta_sd,dt);

  // compute the transition rates
  rate[0] = mu*popsize;		// birth into susceptible class
  rate[1] = (iota+beta*I*dW/dt)/popsize; // force of infection
  rate[2] = mu;			// death from susceptible class
  rate[3] = gamma;	// recovery
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
  incid += trans[3];	// incidence is cumulative recoveries
  if (beta_sd > 0.0) W += (dW-dt)/beta_sd; // increment has mean = 0, variance = dt
"

skel <- "
  double rate[6];		// transition rates
  double term[6];		// transition numbers
  double beta;			// transmission rate
  int k;

  for (k = 0, beta = 0.0; k < nbasis; k++)
     beta += (&beta1)[k]*(&seas_1)[k];

  // compute the transition rates
  rate[0] = mu*popsize;		// birth into susceptible class
  rate[1] = (iota+beta*I)/popsize; // force of infection
  rate[2] = mu;			// death from susceptible class
  rate[3] = gamma;	// recovery
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
"

rinit <- "
  double m = popsize/(S_0+I_0+R_0);
  S = nearbyint(m*S_0);
  I = nearbyint(m*I_0);
  R = nearbyint(m*R_0);
  incid = 0;
  W = 0;
"

library(magrittr)
library(dplyr)

LondonYorke %>%
  filter(town=="New York" & disease=="measles" & year>=1928 & year<=1933) %>%
  select(time,cases) %>%
  pomp(
    times="time",
    t0=1928,
    params=c(
      gamma=26,mu=0.02,iota=0.01,
      beta1=120,beta2=140,beta3=100,
      beta.sd=0.01,
      popsize=5e6,
      rho=0.1,theta=1,
      S.0=0.22,I.0=0.0018,R.0=0.78
    ),
    globals=globals,
    rmeasure=Csnippet(rmeas),
    rprocess=euler.sim(
      step.fun=Csnippet(step.fun),
      delta.t=1/52/20
    ),
    skeleton=vectorfield(Csnippet(skel)),
    covar=covariate_table(
      times=seq(from=1928,to=1934,by=0.01),
      seas=periodic.bspline.basis(
        x=times,
        nbasis=3,
        degree=3,
        period=1
      )
    ),
    partrans=parameter_trans(
      log=c("gamma","mu","iota","beta_sd","theta",sprintf("beta%d",1:3)),
      logit="rho",
      barycentric=c("S_0","I_0","R_0")
    ),
    statenames=c("S","I","R","incid","W"),
    paramnames=c(
      "gamma","mu","iota","beta1","beta2","beta3","beta.sd",
      "popsize","rho","theta","S.0","I.0","R.0"
    ),
    zeronames=c("incid","W"),
    rinit=Csnippet(rinit)
  ) -> po

## compute a trajectory of the deterministic skeleton
X <- trajectory(po,hmax=1/52,format="data.frame")
plot(incid~time,data=X,type='l')

## simulate from the model
x <- simulate(po,nsim=3,format="data.frame")
points(incid~time,data=x,col=as.factor(x$.id),pch=16)

## Compute the likelihood by running a particle filter.
## Note, we must add a 'dmeasure'.
po %>%
  pfilter(
    Np=1000,
    dmeasure=Csnippet("
      lik = dnbinom_mu(cases,theta,rho*incid,give_log);"
    ),
    statenames=c("S","I","R","incid","W"),
    paramnames=c(
      "gamma","mu","iota","beta1","beta.sd",
      "popsize","rho","theta","S.0","I.0","R.0"
    )
  ) %>% logLik()
