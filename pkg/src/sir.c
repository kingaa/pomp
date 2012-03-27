// dear emacs, please treat this as -*- C++ -*-

#include <Rmath.h>

#include "pomp.h"

static double expit (double x) {
  return 1.0/(1.0 + exp(-x));
}

static double logit (double x) {
  return log(x/(1-x));
}

// static double term_time (double t, double b0, double b1) 
// {
//   static double correction = 0.4958904;
//   double day = 365.0 * (t - floor(t));
//   double interm;
//   interm = ((day >= 7.0 && day < 100.0)
// 	    || (day >= 116.0 && day < 200.0) 
// 	    || (day >= 252.0 && day < 300.0) 
// 	    || (day >= 308.0 && day < 356.0)) ? 1.0 : -1.0;
//   return b0*(1.0+b1*interm)/(1.0+b1*correction);
// }

#define GAMMA       (p[parindex[0]]) // recovery rate
#define MU          (p[parindex[1]]) // baseline birth and death rate
#define IOTA        (p[parindex[2]]) // import rate
#define LOGBETA     (p[parindex[3]]) // transmission rate
#define BETA_SD     (p[parindex[4]]) // environmental stochasticity SD in transmission rate
#define POPSIZE     (p[parindex[5]]) // population size
#define RHO         (p[parindex[6]]) // reporting probability
#define NBASIS      (p[parindex[7]]) // number of periodic B-spline basis functions
#define DEGREE      (p[parindex[8]]) // degree of periodic B-spline basis functions
#define PERIOD      (p[parindex[9]]) // period of B-spline basis functions
#define S0         (p[parindex[10]]) // initial fraction of S
#define I0         (p[parindex[11]]) // initial fraction of I
#define R0         (p[parindex[12]]) // initial fraction of R

#define SUSC      (x[stateindex[0]]) // number of susceptibles
#define INFD      (x[stateindex[1]]) // number of infectives
#define RCVD      (x[stateindex[2]]) // number of recovereds
#define CASE      (x[stateindex[3]]) // number of cases (accumulated per reporting period)
#define W         (x[stateindex[4]]) // integrated white noise

#define DSDT    (f[stateindex[0]])
#define DIDT    (f[stateindex[1]])
#define DRDT    (f[stateindex[2]])
#define DCDT    (f[stateindex[3]])
#define DWDT    (f[stateindex[4]])

#define REPORTS   (y[obsindex[0]])

void _sir_par_untrans (double *pt, double *p, int *parindex) 
{
  pt[parindex[0]] = log(GAMMA);
  pt[parindex[1]] = log(MU);
  pt[parindex[2]] = log(IOTA);
  pt[parindex[4]] = log(BETA_SD);
  pt[parindex[6]] = logit(RHO);
  pt[parindex[10]] = log(S0);
  pt[parindex[11]] = log(I0);
  pt[parindex[12]] = log(R0);
}
 
void _sir_par_trans (double *pt, double *p, int *parindex) 
{
  pt[parindex[0]] = exp(GAMMA);
  pt[parindex[1]] = exp(MU);
  pt[parindex[2]] = exp(IOTA);
  pt[parindex[4]] = exp(BETA_SD);
  pt[parindex[6]] = expit(RHO);
  pt[parindex[10]] = exp(S0);
  pt[parindex[11]] = exp(I0);
  pt[parindex[12]] = exp(R0);
}

void _sir_binom_dmeasure (double *lik, double *y, double *x, double *p, int give_log,
			  int *obsindex, int *stateindex, int *parindex, int *covindex,
			  int ncovars, double *covars, double t) {
  double mean, sd;
  double f;
  mean = CASE*RHO;
  sd = sqrt(CASE*RHO*(1-RHO));
  if (REPORTS > 0) {
    f = pnorm(REPORTS+0.5,mean,sd,1,0)-pnorm(REPORTS-0.5,mean,sd,1,0);
  } else {
    f = pnorm(REPORTS+0.5,mean,sd,1,0);
  }
  *lik = (give_log) ? log(f) : f;
}

void _sir_binom_rmeasure (double *y, double *x, double *p, 
			  int *obsindex, int *stateindex, int *parindex, int *covindex,
			  int ncovars, double *covars, double t) {
  double mean, sd;
  double rep;
  mean = CASE*RHO;
  sd = sqrt(CASE*RHO*(1-RHO));
  rep = nearbyint(rnorm(mean,sd));
  REPORTS = (rep > 0) ? rep : 0;
}

// SIR model with Euler multinomial step
// forced transmission (basis functions passed as covariates)
// constant population size as a parameter
// environmental stochasticity on transmission
void _sir_euler_simulator (double *x, const double *p, 
			   const int *stateindex, const int *parindex, const int *covindex,
			   int covdim, const double *covar, 
			   double t, double dt)
{
  int nrate = 6;
  double rate[nrate];		// transition rates
  double trans[nrate];		// transition numbers
  double beta_var;
  double beta;
  double dW;
  int nseas = (int) NBASIS;	// number of seasonal basis functions
  int deg = (int) DEGREE;	// degree of seasonal basis functions
  double seasonality[nseas];

  // untransform the parameters
  beta_var = BETA_SD*BETA_SD;

  if (nseas <= 0) return;
  periodic_bspline_basis_eval(t,PERIOD,deg,nseas,&seasonality[0]);
  beta = exp(dot_product(nseas,&seasonality[0],&LOGBETA));

  //  test to make sure the parameters and state variable values are sane
  if (!(R_FINITE(beta)) || 
      !(R_FINITE(GAMMA)) ||
      !(R_FINITE(MU)) ||
      !(R_FINITE(BETA_SD)) ||
      !(R_FINITE(IOTA)) ||
      !(R_FINITE(POPSIZE)) ||
      !(R_FINITE(SUSC)) ||
      !(R_FINITE(INFD)) ||
      !(R_FINITE(RCVD)) ||
      !(R_FINITE(CASE)) ||
      !(R_FINITE(W)))
    return;

  if (BETA_SD > 0.0) {		// environmental noise is ON
    dW = rgamma(dt/beta_var,beta_var); // gamma noise, mean=dt, variance=(beta_sd^2 dt)
    if (!(R_FINITE(dW))) return;
  } else {			// environmental noise is OFF
    dW = dt;
  }

  // compute the transition rates
  rate[0] = MU*POPSIZE;		// birth into susceptible class
  rate[1] = (IOTA+beta*INFD*dW/dt)/POPSIZE; // force of infection
  rate[2] = MU;			// death from susceptible class
  rate[3] = GAMMA;		// recovery
  rate[4] = MU;			// death from infectious class
  rate[5] = MU; 		// death from recovered class

  // compute the transition numbers
  trans[0] = rpois(rate[0]*dt);	// births are Poisson
  reulermultinom(2,SUSC,&rate[1],dt,&trans[1]);
  reulermultinom(2,INFD,&rate[3],dt,&trans[3]);
  reulermultinom(1,RCVD,&rate[5],dt,&trans[5]);

  // balance the equations
  SUSC += trans[0]-trans[1]-trans[2];
  INFD += trans[1]-trans[3]-trans[4];
  RCVD += trans[3]-trans[5];
  CASE += trans[3];		// cases are cumulative recoveries
  if (BETA_SD > 0.0) {
    W += (dW-dt)/BETA_SD;	// mean zero, variance = dt
  }

}

void _sir_ODE (double *f, double *x, const double *p, 
	       const int *stateindex, const int *parindex, const int *covindex,
	       int covdim, const double *covar, double t) 
{
  int nrate = 6;
  double rate[nrate];		// transition rates
  double term[nrate];		// terms in the equations
  double beta;
  int nseas = (int) NBASIS;	// number of seasonal basis functions
  int deg = (int) DEGREE;	// degree of seasonal basis functions
  double seasonality[nseas];
  
  if (nseas <= 0) return;
  periodic_bspline_basis_eval(t,PERIOD,deg,nseas,&seasonality[0]);
  beta = exp(dot_product(nseas,&seasonality[0],&LOGBETA));

  // compute the transition rates
  rate[0] = MU*POPSIZE;		// birth into susceptible class
  rate[1] = (IOTA+beta*INFD)/POPSIZE; // force of infection
  rate[2] = MU;			// death from susceptible class
  rate[3] = GAMMA;		// recovery
  rate[4] = MU;			// death from infectious class
  rate[5] = MU; 		// death from recovered class

  // compute the several terms
  term[0] = rate[0];
  term[1] = rate[1]*SUSC;
  term[2] = rate[2]*SUSC;
  term[3] = rate[3]*INFD;
  term[4] = rate[4]*INFD;
  term[5] = rate[5]*RCVD;

  // balance the equations
  DSDT = term[0]-term[1]-term[2];
  DIDT = term[1]-term[3]-term[4];
  DRDT = term[3]-term[5];
  DCDT = term[3];		// accumulate the new I->R transitions
  DWDT = 0;			// no noise, so no noise accumulation

}

#undef SUSC
#undef INFD
#undef RCVD
#undef CASE
#undef W

#define SUSC      (x[stateindex[0]]) // number of susceptibles
#define INFD      (x[stateindex[1]]) // number of infectives
#define RCVD      (x[stateindex[2]]) // number of recovereds
#define POPN      (x[stateindex[3]]) // population size
#define CASE      (x[stateindex[4]]) // number of cases (accumulated per reporting period)

double _sir_rates (int j, double t, double *x, double *p,
		   int *stateindex, int *parindex, int *covindex,
		   int ncovar, double *covar) {
  double beta;
  double rate = 0.0;
  int nseas = (int) NBASIS;	// number of seasonal basis functions
  int deg = (int) DEGREE;	// degree of seasonal basis functions
  double seasonality[nseas];

  switch (j) {
  case 1: 			// birth
    rate = MU*POPN;
    break;
  case 2:			// susceptible death
    rate = MU*SUSC;
    break;
  case 3:			// infection
    periodic_bspline_basis_eval(t,PERIOD,deg,nseas,&seasonality[0]);
    beta = exp(dot_product(nseas,&seasonality[0],&LOGBETA));
    rate = (beta*INFD+IOTA)*SUSC/POPSIZE;
    break;
  case 4:			// infected death
    rate = MU*INFD;
    break;
  case 5:			// recovery
    rate = GAMMA*INFD;
    break;
  case 6:			// recovered death
    rate = MU*RCVD;
    break;
  default:
    error("unrecognized rate code %d",j);
    break;
  }
  return rate;
}

