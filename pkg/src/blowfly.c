// dear emacs, please treat this as -*- C++ -*-

#include <Rmath.h>

#include "pomp.h"

#define P       (p[parindex[0]]) // growth rate
#define NZERO   (p[parindex[1]]) // density-dependence parameter
#define DELTA   (p[parindex[2]]) // survival parameter
#define SIGMAP  (p[parindex[3]]) // recruitment noise SD
#define SIGMAD  (p[parindex[4]]) // survivorship noise SD
#define TAU     (p[parindex[5]]) // delay

#define N      (&x[stateindex[0]]) // total population
#define R       (x[stateindex[1]]) // recruits
#define S       (x[stateindex[2]]) // survivors
#define E       (x[stateindex[3]]) // recruitment noise
#define EPS     (x[stateindex[4]]) // survival noise

// Ricker model with log-normal process noise
void _blowfly_model_simulator (double *x, const double *p, 
			       const int *stateindex, const int *parindex, const int *covindex,
			       int covdim, const double *covar, 
			       double t, double dt)
{
  double e = rgammawn(SIGMAP,dt)/dt;
  double eps = rgammawn(SIGMAD,dt)/dt;
  int tau = nearbyint(TAU);
  int k;
  R = rpois(P*N[tau]*exp(-N[tau]/NZERO)*dt*e);
  S = rbinom(N[0],exp(-DELTA*dt*eps));
  E = e;
  EPS = eps;
  for (k = tau; k > 0; k--) N[k] = N[k-1];
  N[0] = R+S;
}

#undef N
#undef R
#undef S
#undef E
#undef EPS

#undef LOG_P
#undef LOG_NZERO
#undef LOG_DELTA
#undef LOG_SIGMAP
#undef LOG_SIGMAD
#undef TAU

