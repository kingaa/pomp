// dear emacs, please treat this as -*- C++ -*-

#include <Rmath.h>

#include "euler.h"
#include "lookup_table.h"

static double expit (double x) {
  return 1.0/(1.0 + exp(-x));
}

static double logit (double x) {
  return log(x/(1-x));
}

static double term_time (double t, double b0, double b1) 
{
  static double correction = 0.52328767123287671233;
  double day = 365.0 * (t - floor(t));
  int tt = (day >= 7.0 && day <= 100.0) 
    || (day >= 115.0 && day <= 200.0) 
    || (day >= 252.0 && day <= 300.0) 
    || (day >= 308.0 && day <= 356.0);
  return b0 * (1.0 + b1 * (-1.0 + 2.0 * ((double) tt))) / (1.0 + correction * b1);
}

#define LOGGAMMA       (p[parindex[0]]) // recovery rate
#define LOGMU          (p[parindex[1]]) // baseline birth and death rate
#define LOGIOTA        (p[parindex[2]]) // import rate
#define LOGBETA        (p[parindex[3]]) // transmission rate
#define LOGBETA_SD     (p[parindex[4]]) // environmental stochasticity SD in transmission rate
#define LOGPOPSIZE     (p[parindex[5]]) // population size

#define SUSC      (x[stateindex[0]]) // number of susceptibles
#define INFD      (x[stateindex[1]]) // number of infectives
#define RCVD      (x[stateindex[2]]) // number of recovereds
#define CASE      (x[stateindex[3]]) // number of cases (accumulated per reporting period)
#define W         (x[stateindex[4]]) // integrated white noise
#define TRANSTNS  (x[stateindex[5]]) // transition numbers

#define SEASBASIS (covar[covindex[0]]) // first column of seasonality basis in lookup table

// SIR model with Euler multinomial step
// forced transmission (basis functions passed as covariates)
// constant population size as a parameter
// environmental stochasticity on transmission
void sir_euler_simulator (double *x, const double *p, 
			  const int *stateindex, const int *parindex, const int *covindex,
			  int covdim, const double *covar, 
			  double t, double dt)
{
  int nrate = 6;
  double rate[nrate];		// transition rates
  double *trans;		// transition numbers
  double gamma, mu, iota, beta_sd, popsize;
  double beta;
  double dW;			// white noise process
  
  // untransform the parameters
  gamma = exp(LOGGAMMA);
  mu = exp(LOGMU);
  iota = exp(LOGIOTA);
  beta_sd = exp(LOGBETA_SD);
  popsize = exp(LOGPOPSIZE);

  // do table lookup to determine seasonal beta
  beta = exp(dot_product(covdim,&SEASBASIS,&LOGBETA));

  // test to make sure the parameters and state variable values are sane
  if (!(R_FINITE(beta)) || 
      !(R_FINITE(gamma)) ||
      !(R_FINITE(mu)) ||
      !(R_FINITE(beta_sd)) ||
      !(R_FINITE(iota)) ||
      !(R_FINITE(popsize)) ||
      !(R_FINITE(SUSC)) ||
      !(R_FINITE(INFD)) ||
      !(R_FINITE(RCVD)) ||
      !(R_FINITE(CASE)) ||
      !(R_FINITE(W)))
    return;

  if (beta_sd > 0.0) {		// environmental noise is ON
    double beta_var = beta_sd*beta_sd;
    dW = rgamma(dt/beta_var,beta_var); // gamma noise, mean=dt, variance=(beta_sd^2 dt)
    if (!(R_FINITE(dW))) return;
  } else {			// environmental noise is OFF
    dW = dt;
  }

  // compute the transition rates
  rate[0] = mu*popsize;		// birth into susceptible class
  rate[1] = (iota+beta*INFD*dW/dt)/popsize; // force of infection
  rate[2] = mu;			// death from susceptible class
  rate[3] = gamma;		// recovery
  rate[4] = mu;			// death from infectious class
  rate[5] = mu; 		// death from recovered class

  // compute the transition numbers
  trans = &TRANSTNS;
  trans[0] = rpois(rate[0]*dt);	// births are Poisson
  reulermultinom(2,SUSC,&rate[1],dt,&trans[1]);
  reulermultinom(2,INFD,&rate[3],dt,&trans[3]);
  reulermultinom(1,RCVD,&rate[5],dt,&trans[5]);
  trans[6] = dW;

  // balance the equations
  SUSC += trans[0]-trans[1]-trans[2];
  INFD += trans[1]-trans[3]-trans[4];
  RCVD += trans[3]-trans[5];
  CASE += trans[3];		// cases are cumulative recoveries
  if (beta_sd > 0.0) {
    W += (dW-dt)/beta_sd;		// mean zero, variance = dt
  }

}

#undef SUSC   
#undef INFD   
#undef RCVD   
#undef CASE   
#undef W      
#undef TRANSTNS

#define SUSC      (x1[stateindex[0]]) // number of susceptibles
#define INFD      (x1[stateindex[1]]) // number of infectives
#define RCVD      (x1[stateindex[2]]) // number of recovereds
#define CASE      (x1[stateindex[3]]) // number of cases (accumulated per reporting period)
#define W         (x1[stateindex[4]]) // integrated white noise
#define TRANSTNS  (x2[stateindex[5]]) // transition numbers

// SIR model with Euler multinomial step
// forced transmission (basis functions passed as covariates)
// constant population size as a parameter
// environmental stochasticity on transmission
void sir_euler_density (double *f, double *x1, double *x2, const double *p, 
			const int *stateindex, const int *parindex, const int *covindex,
			int covdim, const double *covar, 
			double t, double dt)
{
  int nrate = 6;
  double rate[nrate];		// transition rates
  double *trans;		// transition numbers
  double gamma, mu, iota, beta_sd, popsize;
  double beta;
  double dW;			// white noise process
  
  // untransform the parameters
  gamma = exp(LOGGAMMA);
  mu = exp(LOGMU);
  iota = exp(LOGIOTA);
  beta_sd = exp(LOGBETA_SD);
  popsize = exp(LOGPOPSIZE);

  // do table lookup to determine seasonal beta
  beta = exp(dot_product(covdim,&SEASBASIS,&LOGBETA));

  // test to make sure the parameters and state variable values are sane
  if (!(R_FINITE(beta)) || 
      !(R_FINITE(gamma)) ||
      !(R_FINITE(mu)) ||
      !(R_FINITE(beta_sd)) ||
      !(R_FINITE(iota)) ||
      !(R_FINITE(popsize)) ||
      !(R_FINITE(SUSC)) ||
      !(R_FINITE(INFD)) ||
      !(R_FINITE(RCVD)) ||
      !(R_FINITE(CASE)) ||
      !(R_FINITE(W))) {
    *f = R_NaN;
    return;
  }

  // compute the transition rates
  trans = &TRANSTNS;
  dW = trans[6];
  if (beta_sd > 0.0) {		// environmental noise is ON
    double beta_var = beta_sd*beta_sd;
    *f = dgamma(dW,dt/beta_var,beta_var,1);
  } else {			// environmental noise is OFF
    *f = 0;			// THIS ASSUMES THAT dw = dt !!!
  }
  rate[0] = mu*popsize;		// birth into susceptible class
  rate[1] = (iota+beta*INFD*dW/dt)/popsize; // force of infection
  rate[2] = mu;			// death from susceptible class
  rate[3] = gamma;		// recovery
  rate[4] = mu;			// death from infectious class
  rate[5] = mu; 		// death from recovered class

  // compute the transition numbers
  *f += dpois(trans[0],rate[0]*dt,1); // births are Poisson
  *f += deulermultinom(2,SUSC,&rate[1],dt,&trans[1],1);
  *f += deulermultinom(2,INFD,&rate[3],dt,&trans[3],1);
  *f += deulermultinom(1,RCVD,&rate[5],dt,&trans[5],1);
}

#undef GAMMA
#undef MU   
#undef IOTA 
#undef BETA 
#undef BETA_SD
#undef POPSIZE

#undef SUSC   
#undef INFD   
#undef RCVD   
#undef CASE   
#undef W      
#undef TRANSTNS

#undef SEASBASIS
