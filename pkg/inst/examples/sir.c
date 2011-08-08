// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>

static double expit (double x) {
  return 1.0/(1.0 + exp(-x));
}

static double logit (double x) {
  return log(x/(1-x));
}

// some macros to make the code easier to read.
// note how variables and parameters use the indices.
// the integer vectors parindex, stateindex, obsindex
// are computed by pomp by matching the names of input vectors 
// against parnames, statenames, and obsnames, respectively.

#define LOGGAMMA       (p[parindex[0]]) // recovery rate
#define LOGMU          (p[parindex[1]]) // baseline birth and death rate
#define LOGIOTA        (p[parindex[2]]) // import rate
#define LOGBETA        (p[parindex[3]]) // transmission rate
#define LOGBETA_SD     (p[parindex[4]]) // environmental stochasticity SD in transmission rate
#define POPSIZE        (p[parindex[5]]) // population size
#define LOGRHO         (p[parindex[6]]) // reporting probability
#define NBASIS         (p[parindex[7]]) // number of periodic B-spline basis functions
#define DEGREE         (p[parindex[8]]) // degree of periodic B-spline basis functions
#define PERIOD         (p[parindex[9]]) // period of B-spline basis functions

#define SUSC         (x[stateindex[0]]) // number of susceptibles
#define INFD         (x[stateindex[1]]) // number of infectives
#define RCVD         (x[stateindex[2]]) // number of recovereds
#define CASE         (x[stateindex[3]]) // number of cases (accumulated per reporting period)
#define W            (x[stateindex[4]]) // integrated white noise

#define REPORTS        (y[obsindex[0]]) // the data

// the measurement model.
// note that dmeasure and rmeasure are mirrors of each other.

void binomial_dmeasure (double *lik, double *y, double *x, double *p, int give_log,
			  int *obsindex, int *stateindex, int *parindex, int *covindex,
			  int ncovars, double *covars, double t) {
  *lik = dbinom(REPORTS,nearbyint(CASE),exp(LOGRHO),give_log);
}

void binomial_rmeasure (double *y, double *x, double *p, 
			  int *obsindex, int *stateindex, int *parindex, int *covindex,
			  int ncovars, double *covars, double t) {
  REPORTS = rbinom(nearbyint(CASE),exp(LOGRHO));
}

#undef REPORTS

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
  double gamma, mu, iota, beta_sd, beta_var;
  double beta;			// transmission rate
  double dW;			// white noise increment
  int nseas = (int) NBASIS;	// number of seasonal basis functions
  int deg = (int) DEGREE;	// degree of seasonal basis functions
  double seasonality[nseas];
  void (*eval_basis)(double,double,int,int,double*);
  void (*reulmult)(int,double,double*,double,double*);
  double (*dot)(int,const double *, const double *);

  // untransform the parameters
  gamma = exp(LOGGAMMA);
  mu = exp(LOGMU);
  iota = exp(LOGIOTA);
  beta_sd = exp(LOGBETA_SD);
  beta_var = beta_sd*beta_sd;

  // to evaluate the basis functions and compute the transmission rate, use some of 
  // pomp's built-in C-level facilities:
  eval_basis = (void (*)(double,double,int,int,double*)) R_GetCCallable("pomp","periodic_bspline_basis_eval");
  dot = (double (*)(int,const double *,const double *)) R_GetCCallable("pomp","dot_product");
  // pomp's C-level eulermultinomial simulator
  reulmult = (void (*)(int,double,double*,double,double*)) R_GetCCallable("pomp","reulermultinom");

  // compute transmission rate from seasonality
  if (nseas <= 0) return;
  (*eval_basis)(t,PERIOD,deg,nseas,&seasonality[0]); // evaluate the periodic B-spline basis
  beta = exp((*dot)(nseas,&seasonality[0],&LOGBETA)); // compute the transmission rate

  // test to make sure the parameters and state variable values are sane
  // if not, return without doing any work
  if (!(R_FINITE(beta)) || 
      !(R_FINITE(gamma)) ||
      !(R_FINITE(mu)) ||
      !(R_FINITE(beta_sd)) ||
      !(R_FINITE(iota)) ||
      !(R_FINITE(POPSIZE)) ||
      !(R_FINITE(SUSC)) ||
      !(R_FINITE(INFD)) ||
      !(R_FINITE(RCVD)) ||
      !(R_FINITE(CASE)) ||
      !(R_FINITE(W)))
    return;

  // compute the environmental stochasticity
  if (beta_sd > 0.0) {		// environmental noise is ON
    dW = rgamma(dt/beta_var,beta_var); // gamma noise, mean=dt, variance=(beta_sd^2 dt)
    if (!(R_FINITE(dW))) return;
  } else {			// environmental noise is OFF
    dW = dt;
  }

  // compute the transition rates
  rate[0] = mu*POPSIZE;		// birth into susceptible class
  rate[1] = (iota+beta*INFD*dW/dt)/POPSIZE; // force of infection
  rate[2] = mu;			// death from susceptible class
  rate[3] = gamma;		// recovery
  rate[4] = mu;			// death from infectious class
  rate[5] = mu; 		// death from recovered class

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
  if (beta_sd > 0.0) {
    W += (dW-dt)/beta_sd;	// increment mean zero, variance = dt
  }

}

#define DSDT    (f[stateindex[0]])
#define DIDT    (f[stateindex[1]])
#define DRDT    (f[stateindex[2]])
#define DCDT    (f[stateindex[3]])

void sir_ODE (double *f, double *x, const double *p, 
	      const int *stateindex, const int *parindex, const int *covindex,
	      int covdim, const double *covar, double t) 
{
  int nrate = 6;
  double rate[nrate];		// transition rates
  double term[nrate];		// terms in the equations
  double gamma, mu, iota;
  double beta;
  int nseas = (int) NBASIS;	// number of seasonal basis functions
  int deg = (int) DEGREE;	// degree of seasonal basis functions
  double seasonality[nseas];
  void (*eval_basis)(double,double,int,int,double*);
  double (*dot)(int,const double *, const double *);
  
  // untransform the parameters
  gamma = exp(LOGGAMMA);
  mu = exp(LOGMU);
  iota = exp(LOGIOTA);

  // to evaluate the basis functions and compute the transmission rate, use some of 
  // pomp's built-in C-level facilities:
  eval_basis = (void (*)(double,double,int,int,double*)) R_GetCCallable("pomp","periodic_bspline_basis_eval");
  dot = (double (*)(int,const double *,const double *)) R_GetCCallable("pomp","dot_product");

  // compute transmission rate from seasonality
  if (nseas <= 0) return;
  (*eval_basis)(t,PERIOD,deg,nseas,&seasonality[0]); // evaluate the periodic B-spline basis
  beta = exp((*dot)(nseas,&seasonality[0],&LOGBETA));

  // compute the transition rates
  rate[0] = mu*POPSIZE;		// birth into susceptible class
  rate[1] = (iota+beta*INFD)/POPSIZE; // force of infection
  rate[2] = mu;			// death from susceptible class
  rate[3] = gamma;		// recovery
  rate[4] = mu;			// death from infectious class
  rate[5] = mu; 		// death from recovered class

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
  // note that we cannot do the same accumulation trick using 'zeronames' 
  // that we do in the stochastic simulator
  // instead, we'll try to get roughly the number of new infections in a week using:
  DCDT = DIDT*gamma/52;	// this assumes that the reporting rate is weekly!

}

#undef DSDT
#undef DIDT
#undef DRDT
#undef DCDT

#undef SUSC
#undef INFD
#undef RCVD
#undef CASE
#undef W

#undef LOGGAMMA
#undef LOGMU
#undef LOGIOTA
#undef LOGBETA
#undef LOGBETA_SD
#undef POPSIZE
#undef LOGRHO
#undef NBASIS
#undef DEGREE
#undef PERIOD
