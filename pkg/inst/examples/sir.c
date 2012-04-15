// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <pomp.h>

// SIR example as described in the "Advanced Topics in pomp" vignette.
// for a demonstration of how to compile, link, and run this example,
// do 'demo("sir",package="pomp")' at the R prompt.

// some macros to make the code easier to read.
// note how variables and parameters use the indices.
// the integer vectors parindex, stateindex, obsindex
// are computed by pomp by matching the names of input vectors 
// against parnames, statenames, and obsnames, respectively.

#define GAMMA       (p[parindex[0]]) // recovery rate
#define MU          (p[parindex[1]]) // baseline birth and death rate
#define IOTA        (p[parindex[2]]) // import rate
#define BETA       (&p[parindex[3]]) // transmission rate
#define BETA_SD     (p[parindex[4]]) // environmental stochasticity SD in transmission rate
#define POPSIZE     (p[parindex[5]]) // population size
#define RHO         (p[parindex[6]]) // reporting probability
#define NBASIS      (p[parindex[7]]) // number of periodic B-spline basis functions
#define DEGREE      (p[parindex[8]]) // degree of periodic B-spline basis functions
#define PERIOD      (p[parindex[9]]) // period of B-spline basis functions

#define SUSC        (x[stateindex[0]]) // number of susceptibles
#define INFD        (x[stateindex[1]]) // number of infectives
#define RCVD        (x[stateindex[2]]) // number of recovereds
#define CASE        (x[stateindex[3]]) // number of cases (accumulated per reporting period)
#define W           (x[stateindex[4]]) // integrated white noise

#define REPORTS     (y[obsindex[0]]) // the data

#define DSDT        (f[stateindex[0]])
#define DIDT        (f[stateindex[1]])
#define DRDT        (f[stateindex[2]])
#define DCDT        (f[stateindex[3]])
#define DWDT        (f[stateindex[4]])

// the measurement model.
// note that dmeasure and rmeasure are mirrors of each other.

void binomial_dmeasure (double *lik, double *y, double *x, double *p, int give_log,
			  int *obsindex, int *stateindex, int *parindex, int *covindex,
			  int ncovars, double *covars, double t) {
  *lik = dbinom(REPORTS,nearbyint(CASE),RHO,give_log);
}

void binomial_rmeasure (double *y, double *x, double *p, 
			  int *obsindex, int *stateindex, int *parindex, int *covindex,
			  int ncovars, double *covars, double t) {
  REPORTS = rbinom(nearbyint(CASE),RHO);
}

// the process model:
// an SIR model with Euler-multinomial step,
// transmission is seasonal, implemented using B-spline basis functions passed as covariates.
// the population size is constant and specified by a parameter.
// there is environmental stochasticity in transmission (dW).
void sir_euler_simulator (double *x, const double *p, 
			  const int *stateindex, const int *parindex, const int *covindex,
			  int covdim, const double *covar, 
			  double t, double dt)
{
  int nrate = 6;
  double rate[nrate];		// transition rates
  double trans[nrate];		// transition numbers
  double beta;			// transmission rate
  double dW;			// white noise increment
  int nbasis = (int) NBASIS;	// number of seasonal basis functions
  int deg = (int) DEGREE;	// degree of seasonal basis functions
  int k;
  double seasonality[nbasis];
  void (*eval_basis)(double,double,int,int,double*);
  void (*reulmult)(int,double,double*,double,double*);

  // to evaluate the basis functions and compute the transmission rate, use some of 
  // pomp's built-in C-level facilities:
  eval_basis = (void (*)(double,double,int,int,double*)) R_GetCCallable("pomp","periodic_bspline_basis_eval");
  // pomp's C-level eulermultinomial simulator
  reulmult = (void (*)(int,double,double*,double,double*)) R_GetCCallable("pomp","reulermultinom");

  // compute transmission rate from seasonality
  if (nbasis <= 0) return;
  eval_basis(t,PERIOD,deg,nbasis,&seasonality[0]); // evaluate the periodic B-spline basis
  for (k = 0, beta = 0; k < nbasis; k++) 
    beta += BETA[k]*seasonality[k];

  // compute the environmental stochasticity
  dW = rgammawn(BETA_SD,dt);

  // compute the transition rates
  rate[0] = MU*POPSIZE;		// birth into susceptible class
  rate[1] = (IOTA+beta*INFD*dW/dt)/POPSIZE; // force of infection
  rate[2] = MU;			// death from susceptible class
  rate[3] = GAMMA;		// recovery
  rate[4] = MU;			// death from infectious class
  rate[5] = MU; 		// death from recovered class

  // compute the transition numbers
  trans[0] = rpois(rate[0]*dt);	// births are Poisson
  (*reulmult)(2,SUSC,&rate[1],dt,&trans[1]); // euler-multinomial exits from S class
  (*reulmult)(2,INFD,&rate[3],dt,&trans[3]); // euler-multinomial exits from I class
  (*reulmult)(1,RCVD,&rate[5],dt,&trans[5]); // euler-multinomial exits from R class

  // balance the equations
  SUSC += trans[0]-trans[1]-trans[2];
  INFD += trans[1]-trans[3]-trans[4];
  RCVD += trans[3]-trans[5];
  CASE += trans[3];		// cases are cumulative recoveries
  if (BETA_SD > 0.0) W += (dW-dt)/BETA_SD; // increment has mean = 0, variance = dt

}

void sir_ODE (double *f, double *x, const double *p, 
	      const int *stateindex, const int *parindex, const int *covindex,
	      int covdim, const double *covar, double t) 
{
  int nrate = 6;
  double rate[nrate];		// transition rates
  double term[nrate];		// terms in the equations
  double beta;
  int nbasis = (int) NBASIS;	// number of seasonal basis functions
  int deg = (int) DEGREE;	// degree of seasonal basis functions
  int k;
  double seasonality[nbasis];
  void (*eval_basis)(double,double,int,int,double*);
  
  // to evaluate the basis functions and compute the transmission rate, use some of 
  // pomp's built-in C-level facilities:
  eval_basis = (void (*)(double,double,int,int,double*)) R_GetCCallable("pomp","periodic_bspline_basis_eval");

  // compute transmission rate from seasonality
  if (nbasis <= 0) return;
  eval_basis(t,PERIOD,deg,nbasis,&seasonality[0]); // evaluate the periodic B-spline basis
  for (k = 0, beta = 0; k < nbasis; k++) 
    beta += BETA[k]*seasonality[k];

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

  // assemble the differential equations
  DSDT = term[0]-term[1]-term[2];
  DIDT = term[1]-term[3]-term[4];
  DRDT = term[3]-term[5];
  DCDT = term[3];		// accumulate the new I->R transitions
  DWDT = 0;

}
