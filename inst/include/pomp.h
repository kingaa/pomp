// -*- C++ -*-

#ifndef _POMP_H_
#define _POMP_H_

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

// Prototype for one-step Euler simulator, as used by "euler.simulate":
typedef void euler_step_sim(double *x, const double *p, 
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
// ncovars    = number of covariates
// covars     = pointer to a vector containing the values of the covariates at time t, as interpolated 
//                from the covariate table supplied to 'euler.simulator'
// t          = time at the beginning of the Euler step
// dt         = size (duration) of the Euler step
//  on output:
// x          = contains the new state vector (i.e., at time t+dt)

// Prototype for one-step Euler PDF, as used by "euler.density":
typedef void euler_step_pdf(double *f, 
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

// prototypes for C-level access to Euler-multinomial distribution functions
void reulermultinom (int ntrans, double size, double *rate, double dt, double *trans);
double deulermultinom (int ntrans, double size, double *rate, double dt, double *trans, int give_log);


// Prototype for deterministic skeleton evaluation
typedef void pomp_vectorfield_map (double *f, double *x, double *p, 
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
					   int *stateindex, int *parindex, int *covindex, int *obsindex, 
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
// obsindex   = pointer to vector of integers indexing the variables in 'data' in the order specified by 
//                the 'obsnames' slot
// ncovars    = number of covariates
// covars     = pointer to a vector containing the values of the covariates at time t, as interpolated 
//                from the covariate table supplied to 'pomp.skeleton'
// t          = time at the beginning of the Euler step
//  on output:
// y          = pointer to vector containing simulated observations (length = nobs = nrow(data))


// Prototype for measurement model density evaluator
typedef void pomp_measure_model_density (double *lik, double *y, double *x, double *p, int give_log,
					 int *stateindex, int *parindex, 
					 int *covindex, int *obsindex,
					 int ncovars, double *covars, double t);
// Description:
//  on input:
// y          = pointer to vector of observables at time t
// x          = pointer to state vector at time t
// p          = pointer to parameter vector
// give_log   = should the log likelihood be returned?
// stateindex = pointer to vector of integers indexing the states in 'x' in the order specified by 
//                the 'statenames' slot
// parindex   = pointer to vector of integers indexing the parameters in 'p' in the order specified by 
//                the 'paramnames' slot
// covindex   = pointer to vector of integers indexing the parameters in 'covar'' in the order specified by 
//                the 'covarnames' slot
// obsindex   = pointer to vector of integers indexing the variables in 'data' in the order specified by 
//                the 'obsnames' slot
// ncovars    = number of covariates
// covars     = pointer to a vector containing the values of the covariates at time t, as interpolated 
//                from the covariate table supplied to 'pomp.skeleton'
// t          = time at the beginning of the Euler step
//  on output:
// lik        = pointer to scalar containing (log) likelihood


// lookup-table structure, as used internally
struct lookup_table {
  int length, width;
  int index;
  double *x;
  double *y;
};

// simple linear interpolation of the lookup table (with derivative if desired)
// setting dydt = 0 in the call to 'table_lookup' will bypass computation of the derivative
void table_lookup (struct lookup_table *tab, double x, double *y, double *dydt);

// facility for dotting a vector of parameters ('coef') against a vector of basis-function values ('basis')
double dot_product (int dim, const double *basis, const double *coef);

#endif
