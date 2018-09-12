// -*- C++ -*-

#ifndef _POMP_H_
#define _POMP_H_

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

typedef void periodic_bspline_basis_eval_t (double x, double period, int degree, int nbasis, double *y);
typedef void periodic_bspline_basis_eval_deriv_t (double x, double period, int degree, int nbasis, int deriv, double *y);
typedef const SEXP get_userdata_t (const char *name);
typedef const int *get_userdata_int_t (const char *name);
typedef const double *get_userdata_double_t (const char *name);

// UTILITY FOR GAMMA WHITENOISE
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

// VECTOR DOT-PRODUCT UTILITY
// facility for computing the inner product of
// a vector of parameters ('coef') against a vector of basis-function values ('basis')
static R_INLINE double dot_product (int dim, const double *basis, const double *coef) {
  int j;
  double trans = 0.0;
  for (j = 0; j < dim; j++)
    trans += coef[j]*basis[j];
  return(trans);
}


// LOGIT AND INVERSE LOGIT TRANSFORMATIONS
static R_INLINE double logit (double p) {
  return log(p/(1.0-p));
}

static R_INLINE double expit (double x) {
  return 1.0/(1.0+exp(-x));
}

// C-LEVEL DEFINITIONS OF EULER-MULTINOMIAL DISTRIBUTION FUNCTIONS

// reulermultinom: simulate Euler-multinomial transitions
// Description:
//  on input:
// m      = (positive integer) number of potential transitions ("deaths")
// size   = (positive integer) number of individuals at risk
// rate   = pointer to vector of transition ("death") rates
// dt     = (positive real) duration of time interval
//  on output:
// trans  = pointer to vector containing the random deviates
//          (numbers of individuals making the respective transitions)
// See '?reulermultinom' and vignettes for more on the Euler-multinomial
// distributions.
//
// NB: 'reulermultinom' does not call GetRNGstate() and PutRNGstate() internally
// this must be done by the calling program
// But note that when reulermultinom is called inside a pomp 'rprocess',
// there is no need to call {Get,Put}RNGState() as this is handled by pomp

static R_INLINE void reulermultinom (int m, double size, const double *rate,
  double dt, double *trans) {
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
      warningcall(R_NilValue,"in 'reulermultinom': result of binomial draw is not finite");
    m -= 1;
    for (k = 0; k < m; k++) {
      if (rate[k] > p) p = rate[k];
      trans[k] = ((size > 0) && (p > 0)) ? rbinom(size,rate[k]/p) : 0;
      if (!(R_FINITE(size)&&R_FINITE(p)&&R_FINITE(rate[k])&&R_FINITE(trans[k])))
        warningcall(R_NilValue,"in 'reulermultinom': result of binomial draw is not finite");
      size -= trans[k];
      p -= rate[k];
    }
    trans[m] = size;
  } else {
    for (k = 0; k < m; k++) trans[k] = 0.0;
  }
}

// deulermultinom: probabilities of Euler-multinomial transitions
// Description:
//  on input:
// m        = (positive integer) number of potential transitions ("deaths")
// size     = (positive integer) number of individuals at risk
// rate     = pointer to vector of transition ("death") rates
// dt       = (positive real) duration of time interval
// trans    = pointer to vector containing the data
//            (numbers of individuals making the respective transitions)
// give_log = 1 if log probability is desired; 0 if probability is desired
//  return value: probability or log probability (as requested)
//
// See '?deulermultinom' and vignettes for more on the Euler-multinomial
// distributions.

static R_INLINE double deulermultinom (int m, double size, const double *rate,
  double dt, double *trans, int give_log) {
  double p = 0.0;
  double n = 0.0;
  double ff = 0.0;
  int k;
  if ((dt < 0.0) || (size < 0.0) || (floor(size+0.5) != size)) {
    warningcall(R_NilValue,"in 'deulermultinom': NaNs produced");
    return R_NaN;
  }
  for (k = 0; k < m; k++) {
    if (rate[k] < 0.0) {
      warningcall(R_NilValue,"in 'deulermultinom': NaNs produced");
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
    }
    n -= trans[k];
    p -= rate[k];
  }
  ff = (give_log) ? ff : exp(ff);
  return ff;
}

// dmultinom: multinomial probabilities
// Description:
//  on input:
// m        = (positive integer) dimension of the random variable
// prob     = pointer to m-vector of probabilities
// x        = pointer to m-vector containing the data
// give_log = 1 if log probability is desired; 0 if probability is desired
//  return value: probability or log probability (as requested)

static R_INLINE double dmultinom (int m, const double *prob, double *x, int give_log) {
  double p = 0.0;
  double n = 0.0;
  double ff = 0.0;
  int k;

  for (k = 0; k < m; k++) {
    if (prob[k] < 0.0) {
      warningcall(R_NilValue,"in 'dmultinom': NaNs produced");
      return R_NaN;
    }

    if ((x[k] < 0.0) || (floor(x[k]+0.5) != x[k])) {
      ff = (give_log) ? R_NegInf: 0.0;
      return ff;
    }

    p += prob[k]; // sum of probabilities
    n += x[k]; // total number of events
  }

  for (k = 0; k < m; k++) {
    if ((n > 0) && (p > 0)) {
      if (prob[k] > p) p = prob[k];
      ff += dbinom(x[k],n,prob[k]/p,1);
    } else if (x[k] < 0.0) {
      ff = R_NegInf;
      return ff;
    }

    n -= x[k];
    p -= prob[k];
  }

  ff = (give_log) ? ff : exp(ff);
  return ff;
}

// C-LEVEL DEFINITIONS OF LOG-BARYCENTRIC TRANSFORMATION.
// USEFUL FOR WORKING WITH PARAMETERS CONSTRAINED TO SUM TO 1

// TRANSFORMS TO LOG BARYCENTRIC COORDINATES
// on input:
// x =  pointer to vector of parameters to be tranformed to
//      log barycentric coordinates,
// n =  length of vector.
// on output:
// xt = pointer to vector of log barycentric coordinates
static R_INLINE void to_log_barycentric (double *xt, const double *x, int n) {
  double sum;
  int i;
  for (i = 0, sum = 0.0; i < n; i++) sum += x[i];
  for (i = 0; i < n; i++) xt[i] = log(x[i]/sum);
}

// TRANSFORMS FROM LOG BARYCENTRIC COORDINATES
// on input:
// x =  pointer to vector of parameters in log barycentric coordinates,
// n =  length of vector.
// on output:
// xt = pointer to vector of coordinates on unit simplex
static R_INLINE void from_log_barycentric (double *xt, const double *x, int n) {
  double sum;
  int i;
  for (i = 0, sum = 0.0; i < n; i++) sum += (xt[i] = exp(x[i]));
  for (i = 0; i < n; i++) xt[i] /= sum;
}

// UTILITY FOR EXPONENTIAL/GEOMETRIC RATE CONVERSION
// This function computes r such that if
// N ~ geometric(prob=1-exp(-r dt)) and T ~ exponential(rate=R),
// then E[N dt] = E[T]
// i.e., the rate r for an Euler process that gives the same
// expected waiting time as the exponential process it approximates.
// In particular r -> R as dt -> 0.
static R_INLINE double exp2geom_rate_correction (double R, double dt) {
  return (dt > 0) ? log(1.0+R*dt)/dt : R;
}

static R_INLINE double rbetabinom (double size, double prob, double theta) {
  return rbinom(size,rbeta(prob*theta,(1.0-prob)*theta));
}

static R_INLINE double dbetabinom (double x, double size, double prob,
  double theta, int give_log) {
  double a = theta*prob;
  double b = theta*(1.0-prob);
  double f = lchoose(size,x)-lbeta(a,b)+lbeta(a+x,b+size-x);
  return (give_log) ? f : exp(f);
}

static R_INLINE double rbetanbinom (double mu, double size, double theta) {
  double prob = size/(size+mu);
  return rnbinom(size,rbeta(prob*theta,(1.0-prob)*theta));
}

static R_INLINE double dbetanbinom (double x, double mu, double size,
  double theta, int give_log) {
  double prob = size/(size+mu);
  double a = theta*prob;
  double b = theta*(1.0-prob);
  double f = lchoose(size+x-1,size-1)-lbeta(a,b)+lbeta(a+size,b+x);
  return (give_log) ? f : exp(f);
}

// THE FOLLOWING ARE C PROTOTYPES FOR COMPONENTS OF POMP MODELS
// FOR USE WHEN THE LATTER ARE SPECIFIED USING NATIVE CODES COMPILED INTO
// A MANUALLY-LINKED LIBRARY.
// THEY CANNOT BE USED WITHIN C SNIPPETS.

// PROTOTYPE FOR INITIAL-STATE SAMPLER (rinit)
typedef void pomp_rinit(double *x, const double *p, double t,
  const int *stateindex, const int *parindex, const int *covindex,
  const double *covars);
// Description:
//  on input:
// p          = pointer to parameter vector
// t          = time
// stateindex = pointer to vector of integers pointing to the states in 'x' in the order specified by
//                the 'statenames' argument of 'euler.simulator'
// parindex   = pointer to vector of integers pointing to the parameters in 'p' in the order specified by
//                the 'paramnames' argument of 'euler.simulator'
// covindex   = pointer to vector of integers pointing to the covariates in 'covars' in the order
//                specified by the 'covarnames' argument of 'euler.simulator'
// covars     = pointer to a vector containing the values of the covariates at time t, as interpolated
//                from the covariate table supplied to 'euler.simulator'
//  on output:
// x          = contains the state vector (i.e., at time t)
//
// NB: There is no need to call GetRNGstate() or PutRNGstate() in the body of the user-defined function.
//     The RNG is initialized before any call to this function, and the RNG state is written afterward.
//     Inclusion of these calls in the user-defined function may result in significant slowdown.


// PROTOTYPE FOR STOCHASTIC SIMULATION ALGORITHM REACTION-RATE FUNCTION, AS USED BY "GILLESPIE.SIM":
typedef double pomp_ssa_rate_fn(int event, double t, const double *x, const double *p,
  const int *stateindex, const int *parindex, const int *covindex, const double *covars);
// Description:
//  on input:
// event      = integer specifying the number of the reaction whose rate is desired
//                (the first is event is '1')
// t          = time at which the rates are to be evaluated
// x          = vector of state variables
// p          = vector of parameters
// stateindex = pointer to vector of integers pointing to the states in 'x' in the order specified by
//                the 'statenames' argument of 'SSA.simulator'
// parindex   = pointer to vector of integers pointing to the parameters in 'p' in the order specified by
//                the 'paramnames' argument of 'SSA.simulator'
// covindex   = pointer to vector of integers pointing to the covariates in 'covars' in the order
//                specified by the 'covarnames' argument of 'SSA.simulator'
// covars     = pointer to a vector containing the values of the covariates at time t, as interpolated
//                from the covariate table supplied to 'SSA.simulator'
//  returns the rate of the j-th reaction

// PROTOTYPE FOR ONE-STEP SIMULATOR, AS USED BY "EULER.SIM" AND "ONESTEP.SIM":
typedef void pomp_onestep_sim(double *x, const double *p,
  const int *stateindex, const int *parindex, const int *covindex,
  const double *covars, double t, double dt);
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

// PROTOTYPE FOR ONE-STEP LOG PROBABILITY DENSITY FUNCTION, AS USED BY "ONESTEP.DENS":
typedef void pomp_onestep_pdf(double *loglik,
  const double *x1, const double *x2, double t1, double t2, const double *p,
  const int *stateindex, const int *parindex, const int *covindex,
  const double *covars);
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
// covars     = pointer to a vector containing the values of the covariates at time t, as interpolated
//                from the covariate table supplied to 'euler.density'
//  on output:
// loglik     = pointer to the log probability density (a single scalar)

// PROTOTYPE FOR DETERMINISTIC SKELETON EVALUATION
typedef void pomp_skeleton (double *f, const double *x, const double *p,
  const int *stateindex, const int *parindex, const int *covindex,
  const double *covars, double t);

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
// covars     = pointer to a vector containing the values of the covariates at time t, as interpolated
//                from the covariate table supplied to 'pomp.skeleton'
// t          = time at the beginning of the Euler step
//  on output:
// f          = pointer to value of the map or vectorfield (a vector of the same length as 'x')

// PROTOTYPE FOR MEASUREMENT MODEL SIMULATION
typedef void pomp_measure_model_simulator (double *y, const double *x, const double *p,
  const int *obsindex, const int *stateindex, const int *parindex, const int *covindex,
  const double *covars, double t);
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
// covars     = pointer to a vector containing the values of the covariates at time t, as interpolated
//                from the covariate table supplied to 'pomp.skeleton'
// t          = time at the beginning of the Euler step
//  on output:
// y          = pointer to vector containing simulated observations (length = nobs = nrow(data))
//
// NB: There is no need to call GetRNGstate() or PutRNGstate() in the body of the user-defined function.
//     The RNG is initialized before any call to this function, and the RNG state is written afterward.
//     Inclusion of these calls in the user-defined function may result in significant slowdown.

// PROTOTYPE FOR MEASUREMENT MODEL DENSITY EVALUATOR
typedef void pomp_measure_model_density (double *lik, const double *y, const double *x, const double *p, int give_log,
					 const int *obsindex, const int *stateindex, const int *parindex, const int *covindex,
					 const double *covars, double t);
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
// covars     = pointer to a vector containing the values of the covariates at time t, as interpolated
//                from the covariate table supplied to 'pomp.skeleton'
// t          = time at the beginning of the Euler step
//  on output:
// lik        = pointer to scalar containing (log) likelihood

// PROTOTYPE FOR PRIOR SIMULATION
typedef void pomp_rprior (double *p, const int *parindex);
// Description:
//  on input:
// p          = pointer to parameter vector
// parindex   = pointer to vector of integers indexing the parameters in 'p' in the order specified by
//                the 'paramnames' slot
//  on output:
// p          = pointer to vector containing draws from the prior
//
// NB: There is no need to call GetRNGstate() or PutRNGstate() in the body of the user-defined function.
//     The RNG is initialized before any call to this function, and the RNG state is written afterward.
//     Inclusion of these calls in the user-defined function may result in significant slowdown.

// PROTOTYPE FOR PRIOR DENSITY EVALUATION
typedef void pomp_dprior (double *lik, const double *p, int give_log, const int *parindex);
// Description:
//  on input:
// p          = pointer to parameter vector
// give_log   = should the log likelihood be returned?
// parindex   = pointer to vector of integers indexing the parameters
//              in 'p' in the order specified by the 'paramnames' slot
//  on output:
// lik        = pointer to vector containing likelihoods

// PROTOTYPE FOR PARAMETER TRANSFORMATION FUNCTION.
typedef void pomp_transform_fn (double *pt, const double *p, const int *parindex);
// Description:
//  on input:
// p          = pointer to parameter vector
// parindex   = pointer to vector of integers indexing the parameters
//              in 'p' in the order specified by the 'paramnames' slot
//  on output:
// pt         = pointer to transformed parameter vector

#endif
