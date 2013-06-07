// -*- C++ -*-

#ifndef _POMP_H_
#define _POMP_H_

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// facilities for extracting R objects from the 'userdata' slot
const SEXP get_pomp_userdata (const char *name);
const int *get_pomp_userdata_int (const char *name);
const double *get_pomp_userdata_double (const char *name);

// facility for evaluating a set of periodic bspline basis functions
void periodic_bspline_basis_eval (double x, double period, int degree, int nbasis, double *y);

// Prototype for parameter transformation function.
typedef void pomp_transform_fn (double *pt, double *p, int *parindex);

// Prototype for stochastic simulation algorithm reaction-rate function, as used by "gillespie.sim":
typedef double pomp_ssa_rate_fn(int j, double t, const double *x, const double *p,
				int *stateindex, int *parindex, int *covindex,
				int ncovar, double *covars);
// Description:
//  on input:
// j          = integer specifying the number of the reaction whose rate is desired
// t          = time at which the rates are to be evaluated
// x          = vector of state variables
// p          = vector of parameters
// stateindex = pointer to vector of integers pointing to the states in 'x' in the order specified by 
//                the 'statenames' argument of 'SSA.simulator'
// parindex   = pointer to vector of integers pointing to the parameters in 'p' in the order specified by 
//                the 'paramnames' argument of 'SSA.simulator'
// covindex   = pointer to vector of integers pointing to the covariates in 'covars' in the order 
//                specified by the 'covarnames' argument of 'SSA.simulator'
// ncovars    = number of covariates
// covars     = pointer to a vector containing the values of the covariates at time t, as interpolated 
//                from the covariate table supplied to 'SSA.simulator'
//  returns the rate of the j-th reaction

// Prototype for one-step simulator, as used by "euler.sim" and "onestep.sim":
typedef void pomp_onestep_sim(double *x, const double *p, 
			      const int *stateindex, const int *parindex, const int *covindex,
			      int ncovars, const double *covars,
			      double t, double dt);
// Description:
//  on input:
// x          = pointer to state vector
// p          = pointer to parameter vector
// stateindex = pointer to vector of integers pointing to the states in 'x' in the order specified by 
//                the 'statenames' argument of 'euler.simulator'
// parindex   = pointer to vector of integers pointing to the parameters in 'p' in the order specified by 
//                the 'paramnames' argument of 'euler.simulator'
// covindex   = pointer to vector of integers pointing to the covariates in 'covars' in the order 
//                specified by the 'covarnames' argument of 'euler.simulator'
// ncovars    = number of covariates
// covars     = pointer to a vector containing the values of the covariates at time t, as interpolated 
//                from the covariate table supplied to 'euler.simulator'
// t          = time at the beginning of the Euler step
// dt         = size (duration) of the Euler step
//  on output:
// x          = contains the new state vector (i.e., at time t+dt)
//
// NB: There is no need to call GetRNGstate() or PutRNGstate() in the body of the user-defined function.
//     The RNG is initialized before any call to this function, and the RNG state is written afterward.
//     Inclusion of these calls in the user-defined function may result in significant slowdown.


// Prototype for one-step log probability density function, as used by "onestep.dens":
typedef void pomp_onestep_pdf(double *f, 
			      double *x1, double *x2, double t1, double t2, const double *p, 
			      const int *stateindex, const int *parindex, const int *covindex,
			      int ncovars, const double *covars);
// Description:
//  on input:
// x1         = pointer to state vector at time t1
// x2         = pointer to state vector at time t2
// t1         = time corresponding to x1
// t2         = time corresponding to x2
// p          = pointer to parameter vector
// stateindex = pointer to vector of integers indexing the states in 'x' in the order specified by 
//                the 'statenames' argument of 'euler.density'
// parindex   = pointer to vector of integers indexing the parameters in 'p' in the order specified by 
//                the 'paramnames' argument of 'euler.density'
// covindex   = pointer to vector of integers indexing the parameters in 'covar'' in the order specified by 
//                the 'covarnames' argument of 'euler.density'
// ncovars    = number of covariates
// covars     = pointer to a vector containing the values of the covariates at time t, as interpolated 
//                from the covariate table supplied to 'euler.density'
//  on output:
// f          = pointer to the probability density (a single scalar)

// Prototype for deterministic skeleton evaluation
typedef void pomp_skeleton (double *f, double *x, double *p, 
			    int *stateindex, int *parindex, int *covindex, 
			    int ncovars, double *covars, double t);

// Description:
//  on input:
// x          = pointer to state vector at time t
// p          = pointer to parameter vector
// stateindex = pointer to vector of integers indexing the states in 'x' in the order specified by 
//                the 'statenames' slot
// parindex   = pointer to vector of integers indexing the parameters in 'p' in the order specified by 
//                the 'paramnames' slot
// covindex   = pointer to vector of integers indexing the parameters in 'covar'' in the order specified by 
//                the 'covarnames' slot
// ncovars    = number of covariates
// covars     = pointer to a vector containing the values of the covariates at time t, as interpolated 
//                from the covariate table supplied to 'pomp.skeleton'
// t          = time at the beginning of the Euler step
//  on output:
// f          = pointer to value of the map or vectorfield (a vector of the same length as 'x')

// Prototype for measurement model simulation
typedef void pomp_measure_model_simulator (double *y, double *x, double *p, 
					   int *obsindex, int *stateindex, int *parindex, int *covindex,
					   int ncovars, double *covars, double t);
// Description:
//  on input:
// x          = pointer to state vector at time t
// p          = pointer to parameter vector
// obsindex   = pointer to vector of integers indexing the variables in 'y' in the order specified by 
//                the 'obsnames' slot
// stateindex = pointer to vector of integers indexing the states in 'x' in the order specified by 
//                the 'statenames' slot
// parindex   = pointer to vector of integers indexing the parameters in 'p' in the order specified by 
//                the 'paramnames' slot
// covindex   = pointer to vector of integers indexing the parameters in 'covar'' in the order specified by 
//                the 'covarnames' slot
// ncovars    = number of covariates
// covars     = pointer to a vector containing the values of the covariates at time t, as interpolated 
//                from the covariate table supplied to 'pomp.skeleton'
// t          = time at the beginning of the Euler step
//  on output:
// y          = pointer to vector containing simulated observations (length = nobs = nrow(data))
//
// NB: There is no need to call GetRNGstate() or PutRNGstate() in the body of the user-defined function.
//     The RNG is initialized before any call to this function, and the RNG state is written afterward.
//     Inclusion of these calls in the user-defined function may result in significant slowdown.


// Prototype for measurement model density evaluator
typedef void pomp_measure_model_density (double *lik, double *y, double *x, double *p, int give_log,
					 int *obsindex, int *stateindex, int *parindex, int *covindex,
					 int ncovars, double *covars, double t);
// Description:
//  on input:
// y          = pointer to vector of observables at time t
// x          = pointer to state vector at time t
// p          = pointer to parameter vector
// give_log   = should the log likelihood be returned?
// obsindex   = pointer to vector of integers indexing the variables in 'y' in the order specified by 
//                the 'obsnames' slot
// stateindex = pointer to vector of integers indexing the states in 'x' in the order specified by 
//                the 'statenames' slot
// parindex   = pointer to vector of integers indexing the parameters in 'p' in the order specified by 
//                the 'paramnames' slot
// covindex   = pointer to vector of integers indexing the parameters in 'covar'' in the order specified by 
//                the 'covarnames' slot
// ncovars    = number of covariates
// covars     = pointer to a vector containing the values of the covariates at time t, as interpolated 
//                from the covariate table supplied to 'pomp.skeleton'
// t          = time at the beginning of the Euler step
//  on output:
// lik        = pointer to scalar containing (log) likelihood

// This function computes r such that if
// N ~ geometric(prob=1-exp(-r dt)) and T ~ exponential(rate=R),
// then E[N dt] = E[T]
// i.e., the rate r for an Euler process that gives the same
// expected waiting time as the exponential process it approximates.
// In particular r -> R as dt -> 0.
static R_INLINE double exp2geom_rate_correction (double R, double dt) {
  return (dt > 0) ? log(1.0+R*dt)/dt : R;
}

// This function draws a random increment of a gamma whitenoise process.
// This will have expectation=dt and variance=(sigma^2*dt)
// If dW = rgammawn(sigma,dt), then 
// mu dW/dt is a candidate for a random rate process within an
// Euler-multinomial context, i.e., 
// E[mu*dW] = mu*dt and Var[mu*dW] = mu*sigma^2*dt
static R_INLINE double rgammawn (double sigma, double dt) {
  double sigmasq;
  sigmasq = sigma*sigma;
  return (sigmasq > 0) ? rgamma(dt/sigmasq,sigmasq) : dt;
}

// facility for computing the inner product of 
// a vector of parameters ('coef') against a vector of basis-function values ('basis')
static R_INLINE double dot_product (int dim, const double *basis, const double *coef) {
  int j;
  double trans = 0.0;
  for (j = 0; j < dim; j++)
    trans += coef[j]*basis[j];
  return(trans);
}

static R_INLINE double logit (double p) {
  return log(p/(1.0-p));
}

static R_INLINE double expit (double x) {
  return 1.0/(1.0+exp(-x));
}

// C-level definitions of Euler-multinomial distribution functions

// simulate Euler-multinomial transitions
// NB: 'reulermultinom' does not call GetRNGstate() and PutRNGstate() internally
// this must be done by the calling program
// But note that when reulermultinom is called inside a pomp 'rprocess', there is no need to call
// {Get,Put}RNGState() as this is handled by pomp
static void reulermultinom (int m, double size, double *rate, double dt, double *trans) {
  double p = 0.0;
  int j, k;
  if ((size < 0.0) || (dt < 0.0) || (floor(size+0.5) != size)) {
    for (k = 0; k < m; k++) trans[k] = R_NaN;
    return;
  }  
  for (k = 0; k < m; k++) {
    if (rate[k] < 0.0) {
      for (j = 0; j < m; j++) trans[j] = R_NaN;
      return;
    }
    p += rate[k]; // total event rate
  }
  if (p > 0.0) {
    size = rbinom(size,1-exp(-p*dt)); // total number of events
    if (!(R_FINITE(size)))
      warning("reulermultinom: result of binomial draw is not finite");
    m -= 1;
    for (k = 0; k < m; k++) {
      if (rate[k] > p) p = rate[k];
      trans[k] = ((size > 0) && (p > 0)) ? rbinom(size,rate[k]/p) : 0;
      if (!(R_FINITE(size)&&R_FINITE(p)&&R_FINITE(rate[k])&&R_FINITE(trans[k])))
	warning("reulermultinom: result of binomial draw is not finite");
      size -= trans[k];
      p -= rate[k];
    }
    trans[m] = size;
  } else {
    for (k = 0; k < m; k++) trans[k] = 0.0;
  }
}

// compute probabilities of Euler-multinomial transitions
static double deulermultinom (int m, double size, double *rate, double dt, double *trans, int give_log) {
  double p = 0.0;
  double n = 0.0;
  double ff = 0.0;
  int k;
  if ((dt < 0.0) || (size < 0.0) || (floor(size+0.5) != size)) {
    warning("NaNs produced");
    return R_NaN;
  }
  for (k = 0; k < m; k++) {
    if (rate[k] < 0.0) {
      warning("NaNs produced");
      return R_NaN;
    }
    if (trans[k] < 0.0) {
      ff = (give_log) ? R_NegInf: 0.0;
      return ff;
    }
    p += rate[k]; // total event rate
    n += trans[k]; // total number of events
  }
  if (n > size) {
    ff = (give_log) ? R_NegInf: 0.0;
    return ff;
  }
  ff = dbinom(n,size,1-exp(-p*dt),1); // total number of events
  m -= 1;
  for (k = 0; k < m; k++) {
    if ((n > 0) && (p > 0)) {
      if (rate[k] > p) p = rate[k];
      ff += dbinom(trans[k],n,rate[k]/p,1);
    } else if (trans[k] > 0.0) {
      ff = R_NegInf;
      return ff;
    }
    n -= trans[k];
    p -= rate[k];
  }
  ff = (give_log) ? ff : exp(ff);
  return ff;
}

static R_INLINE double rbetabinom (double size, double prob, double theta) {
  return rbinom(size,rbeta(prob*theta,(1.0-prob)*theta));
}

static R_INLINE double dbetabinom (double x, double size, double prob, double theta, int give_log) {
  double a = theta*prob;
  double b = theta*(1.0-prob);
  double f = lchoose(size,x)-lbeta(a,b)+lbeta(a+x,b+size-x);
  return (give_log) ? f : exp(f);
}

static R_INLINE double rbetanbinom (double mu, double size, double theta) {
  double prob = size/(size+mu);
  return rnbinom(size,rbeta(prob*theta,(1.0-prob)*theta));
}

static R_INLINE double dbetanbinom (double x, double mu, double size, double theta, int give_log) {
  double prob = size/(size+mu);
  double a = theta*prob;
  double b = theta*(1.0-prob);
  double f = lchoose(size+x-1,size-1)-lbeta(a,b)+lbeta(a+size,b+x);
  return (give_log) ? f : exp(f);
}

#endif
