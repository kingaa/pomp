// dear emacs, please treat this as -*- C++ -*-

#include <Rmath.h>

#include "pomp.h"

#define R       (p[parindex[0]]) // growth rate
#define K       (p[parindex[1]]) // carrying capacity
#define SIGMA   (p[parindex[2]]) // process noise level
#define THETA   (p[parindex[3]]) // measurement noise level

#define POP         (y[obsindex[0]])
#define N           (x[stateindex[0]])
#define NPRIME      (f[stateindex[0]])

void _parus_lognormal_dmeasure (double *lik, double *y, double *x, double *p, int give_log,
				int *obsindex, int *stateindex, int *parindex, int *covindex,
				int ncovars, double *covars, double t) {
  *lik = dlnorm(POP,log(N),THETA,give_log);
}

void _parus_lognormal_rmeasure (double *y, double *x, double *p, 
				int *obsindex, int *stateindex, int *parindex, int *covindex,
				int ncovars, double *covars, double t) {
  POP = rlnorm(log(N),THETA);
}

void _parus_poisson_dmeasure (double *lik, double *y, double *x, double *p, int give_log,
			       int *obsindex, int *stateindex, int *parindex, int *covindex,
			       int ncovars, double *covars, double t) {
  *lik = dpois(POP,N+1.0e-10,give_log);
}

void _parus_poisson_rmeasure (double *y, double *x, double *p, 
			      int *obsindex, int *stateindex, int *parindex, int *covindex,
			      int ncovars, double *covars, double t) {
  POP = rpois(N+1.0e-10);
}

void _parus_nbinom_dmeasure (double *lik, double *y, double *x, double *p, int give_log,
			     int *obsindex, int *stateindex, int *parindex, int *covindex,
			     int ncovars, double *covars, double t) {
  *lik = dnbinom_mu(POP,1.0/THETA,N+1.0e-10,give_log);
}

void _parus_nbinom_rmeasure (double *y, double *x, double *p, 
			     int *obsindex, int *stateindex, int *parindex, int *covindex,
			     int ncovars, double *covars, double t) {
  POP = rnbinom_mu(1.0/THETA,N+1.0e-10);
}

void _parus_gompertz_simulator (double *x, const double *p, 
			  const int *stateindex, const int *parindex, const int *covindex,
			  int covdim, const double *covar, 
			  double t, double dt)
{
  double S = exp(-R*dt);
  double eps = (SIGMA > 0.0) ? exp(rnorm(0,SIGMA)) : 1.0;
  N = pow(K,(1-S))*pow(N,S)*eps;
}

// the deterministic skeleton
void _parus_gompertz_skeleton (double *f, double *x, const double *p, 
			       const int *stateindex, const int *parindex, const int *covindex,
			       int covdim, const double *covar, double t) 
{
  double dt = 1.0;
  double S = exp(-R*dt);
  NPRIME = pow(K,(1-S))*pow(N,S);
}

// Ricker model with log-normal process noise
void _parus_ricker_simulator (double *x, const double *p, 
			      const int *stateindex, const int *parindex, const int *covindex,
			      int covdim, const double *covar, 
			      double t, double dt)
{
  double e = (SIGMA > 0.0) ? rnorm(0,SIGMA) : 0.0;
  N = exp(log(N)+R*(1-N/K)+e);
}

void _parus_ricker_skeleton (double *f, double *x, const double *p, 
			     const int *stateindex, const int *parindex, const int *covindex,
			     int covdim, const double *covar, double t) 
{
  NPRIME = exp(log(N)+R*(1-N/K));
}
