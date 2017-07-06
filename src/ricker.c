// dear emacs, please treat this as -*- C++ -*-

#include <Rmath.h>

#include "pomp.h"

#define R       (p[parindex[0]]) // growth rate
#define SIGMA   (p[parindex[1]]) // process noise level
#define PHI     (p[parindex[2]]) // measurement scale parameter
#define C       (p[parindex[3]]) // population scale

#define N       (x[stateindex[0]]) // population size
#define E       (x[stateindex[1]]) // process noise

#define Y       (y[obsindex[0]]) // observed population size

void _ricker_poisson_dmeasure (double *lik, double *y, double *x, double *p, int give_log,
			       int *obsindex, int *stateindex, int *parindex, int *covindex,
			       int ncovars, double *covars, double t) {
  *lik = dpois(Y,PHI*N,give_log);
}

void _ricker_poisson_rmeasure (double *y, double *x, double *p, 
			       int *obsindex, int *stateindex, int *parindex, int *covindex,
			       int ncovars, double *covars, double t) {
  Y = rpois(PHI*N);
}

#undef Y

// Ricker model with log-normal process noise
void _ricker_simulator (double *x, const double *p, 
			const int *stateindex, const int *parindex, const int *covindex,
			int covdim, const double *covar, 
			double t, double dt)
{
  double e = (SIGMA > 0.0) ? rnorm(0,SIGMA) : 0.0;
  N = exp(log(R)+log(N)-C*N+e);
  E = e;
}

void _ricker_skeleton (double *f, double *x, const double *p, 
		       const int *stateindex, const int *parindex, const int *covindex,
		       int covdim, const double *covar, double t) 
{
  f[0] = exp(log(R)+log(N)-C*N);
  f[1] = 0.0;
}

#undef N
#undef E

#undef R
#undef SIGMA
#undef PHI
#undef C
