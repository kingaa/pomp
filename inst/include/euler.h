// -*- C++ -*-

#ifndef _POMP_EULER_H_
#define _POMP_EULER_H_

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

typedef void euler_step_sim(double *x, const double *p, 
			    const int *stateindex, const int *parindex,
			    int ncovar, const double *covar,
			    double t, double dt);
// Description of Euler step functions (type euler_step_sim)
//  on input:
// x          = pointer to state vector
// p          = pointer to parameter vector
// stateindex = pointer to vector of integers pointing to the states in 'x' in the order specified by 
//                the 'statenames' argument of 'euler.simulator'
// parindex   = pointer to vector of integers pointing to the parameters in 'p' in the order specified by 
//                the 'paramnames' argument of 'euler.simulator'
// ncovar     = number of covariates
// covar      = pointer to a vector containing the values of the covariates at time t, as interpolated 
//                from the covariate table supplied to 'euler.simulator'
// t          = time at the beginning of the Euler step
// dt         = size (duration) of the Euler step
//  on output:
// x          = contains the new state vector (i.e., at time t+dt)

typedef void euler_step_pdf(double *f, 
			    double *x1, double *x2, const double *p, 
			    const int *stateindex, const int *parindex,
			    int ncovar, const double *covar,
			    double t, double dt);

void reulermultinom (int ntrans, double size, double *rate, double dt, double *trans);
double deulermultinom (int ntrans, double size, double *rate, double dt, double *trans, int give_log);

#endif
