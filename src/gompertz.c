// dear emacs, please treat this as -*- C++ -*-

#include <Rmath.h>

#include "pomp.h" // in R, do 'system.file("include/pomp.h",package="pomp")' to find this header file

// define some macros to make the code easier to read
#define R       (p[parindex[0]]) // growth rate
#define K       (p[parindex[1]]) // carrying capacity
#define SIGMA   (p[parindex[2]]) // process noise level
#define TAU     (p[parindex[3]]) // measurement noise level

#define Y           (y[obsindex[0]])   // observed population size
#define X           (x[stateindex[0]]) // actual population size
#define XPRIME      (f[stateindex[0]]) // new population size (for skeleton function only)

// normal measurement model density
void _gompertz_normal_dmeasure (double *lik, double *y, double *x, double *p, int give_log,
  int *obsindex, int *stateindex, int *parindex, int *covindex,
  double *covars, double t) {
  *lik = dlnorm(Y,log(X),TAU,give_log);
}

// normal measurement model simulator
void _gompertz_normal_rmeasure (double *y, double *x, double *p,
  int *obsindex, int *stateindex, int *parindex, int *covindex,
  double *covars, double t) {
  Y = rlnorm(log(X),TAU);
}

// stochastic Gompertz model with log-normal process noise
void _gompertz_simulator (double *x, const double *p,
  const int *stateindex, const int *parindex, const int *covindex,
  int covdim, const double *covar, double t, double dt)
{
  double S = exp(-R*dt);
  double eps = (SIGMA > 0.0) ? exp(rnorm(0,SIGMA)) : 1.0;
  X = pow(K,(1-S))*pow(X,S)*eps; // note X is over-written by this line
}

// the deterministic skeleton
void _gompertz_skeleton (double *f, double *x, const double *p,
  const int *stateindex, const int *parindex, const int *covindex,
  const double *covar, double t)
{
  double dt = 1.0;
  double S = exp(-R*dt);
  XPRIME = pow(K,(1-S))*pow(X,S); // X is not over-written in the skeleton function
}

#undef R
#undef K
#undef SIGMA
#undef TAU
#undef Y
#undef X
#undef XPRIME

#define r		(__p[__parindex[0]])
#define K		(__p[__parindex[1]])
#define sigma		(__p[__parindex[2]])
#define tau		(__p[__parindex[3]])
#define X_0		(__p[__parindex[4]])
#define T_r		(__pt[__parindex[0]])
#define T_K		(__pt[__parindex[1]])
#define T_sigma		(__pt[__parindex[2]])
#define T_tau		(__pt[__parindex[3]])
#define T_X_0		(__pt[__parindex[4]])

void _gompertz_to_trans (double *__pt, const double *__p, const int *__parindex)
{
  T_r = log(r);
  T_K = log(K);
  T_sigma = log(sigma);
  T_tau = log(tau);
  T_X_0 = log(X_0);
}

void _gompertz_from_trans (double *__p, const double *__pt, const int *__parindex)
{
  r = exp(T_r);
  K = exp(T_K);
  sigma = exp(T_sigma);
  tau = exp(T_tau);
  X_0 = exp(T_X_0);
}

#undef r
#undef K
#undef sigma
#undef tau
#undef X_0
#undef T_r
#undef T_K
#undef T_sigma
#undef T_tau
#undef T_X_0
