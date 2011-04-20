// dear emacs, please treat this as -*- C++ -*-

#include <Rmath.h>

#include "pomp.h"

#define LOG_R       (p[parindex[0]]) // growth rate
#define LOG_K       (p[parindex[1]]) // carrying capacity
#define LOG_SIGMA   (p[parindex[2]]) // process noise level
#define LOG_TAU     (p[parindex[3]]) // measurement noise level

#define X           (x[stateindex[0]]) // population size

#define Y           (y[obsindex[0]]) // observed population size

void _gompertz_normal_dmeasure (double *lik, double *y, double *x, double *p, int give_log,
				int *obsindex, int *stateindex, int *parindex, int *covindex,
				int ncovars, double *covars, double t) {
  *lik = dlnorm(Y,log(X),exp(LOG_TAU),give_log);
}

void _gompertz_normal_rmeasure (double *y, double *x, double *p, 
				int *obsindex, int *stateindex, int *parindex, int *covindex,
				int ncovars, double *covars, double t) {
  Y = rlnorm(log(X),exp(LOG_TAU));
}

#undef Y

// Ricker model with log-normal process noise
void _gompertz_simulator (double *x, const double *p, 
			  const int *stateindex, const int *parindex, const int *covindex,
			  int covdim, const double *covar, 
			  double t, double dt)
{
  double S = exp(-exp(LOG_R)*dt);
  double sigma = exp(LOG_SIGMA);
  double logeps = (sigma > 0.0) ? rnorm(0,sigma) : 0.0;
  X = exp((1-S)*LOG_K)*pow(X,S)*exp(logeps);
}

void _gompertz_skeleton (double *f, double *x, const double *p, 
			 const int *stateindex, const int *parindex, const int *covindex,
			 int covdim, const double *covar, double t) 
{
  double dt = 1.0;
  double S = exp(-exp(LOG_R)*dt);
  f[0] = exp((1-S)*LOG_K)*pow(X,S);
}

#undef N
#undef E

#undef LOG_K
#undef LOG_R
#undef LOG_SIGMA
#undef LOG_TAU
