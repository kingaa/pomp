// -*- C++ -*-

#ifndef _POMP_H_
#define _POMP_H_

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

// prototypes for C-level access to Euler-multinomial distribution functions
// NB: 'reulermultinom' does not call GetRNGstate() and PutRNGstate() internally
void reulermultinom (int ntrans, double size, double *rate, double dt, double *trans);
double deulermultinom (int ntrans, double size, double *rate, double dt, double *trans, int give_log);

// This function computes r such that if
// N ~ geometric(prob=1-exp(-r dt)) and T ~ exponential(rate=R),
// then E[N dt] = E[T]
// i.e., the rate r for an Euler process that gives the same
// expected waiting time as the exponential process it approximates.
// In particular r -> R as dt -> 0.
inline double exp2geom_rate_correction (double R, double dt) {
  return (dt > 0) ? log(1.0+R*dt)/dt : R;
}

// This function draws a random increment of a gamma whitenoise process.
// This will have expectation=dt and variance=(sigma^2*dt)
// If dW = rgammawn(sigma,dt), then 
// mu dW/dt is a candidate for a random rate process within an
// Euler-multinomial context, i.e., 
// E[mu*dW] = mu*dt and Var[mu*dW] = mu*sigma^2*dt
double rgammawn (double sigma, double dt);

// facility for computing the inner produce of 
// a vector of parameters ('coef') against a vector of basis-function values ('basis')
double dot_product (int dim, const double *basis, const double *coef);

// facility for computing evaluating a basis of periodic bsplines
void periodic_bspline_basis_eval (double x, double period, int degree, int nbasis, double *y);

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

#endif
