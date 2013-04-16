// dear emacs, please treat this as -*- C++ -*-

#include <Rmath.h>

#include "pomp.h"

#define P       (p[parindex[0]]) // growth rate
#define NZERO   (p[parindex[1]]) // density-dependence parameter
#define DELTA   (p[parindex[2]]) // survival parameter
#define SIGMAP  (p[parindex[3]]) // recruitment noise SD
#define SIGMAD  (p[parindex[4]]) // survivorship noise SD
#define SIGMAY  (p[parindex[5]]) // survivorship noise SD

#define N      (&x[stateindex[0]]) // total population
#define R       (x[stateindex[1]]) // recruits
#define S       (x[stateindex[2]]) // survivors
#define E       (x[stateindex[3]]) // recruitment noise
#define EPS     (x[stateindex[4]]) // survival noise

#define Y       (y[obsindex[0]]) // observation

// Ricker model with log-normal process noise
static void _blowfly_simulator (double *x, const double *p, 
				const int *stateindex, const int *parindex, const int *covindex,
				int covdim, const double *covar, 
				double t, double dt, int tau)
{
  double e = rgammawn(SIGMAP,dt)/dt;
  double eps = rgammawn(SIGMAD,dt)/dt;
  int k;

  R = rpois(P*N[tau]*exp(-N[tau]/NZERO)*dt*e);
  S = rbinom(N[0],exp(-DELTA*dt*eps));
  E = e;
  EPS = eps;
  for (k = tau; k > 0; k--) N[k] = N[k-1];
  N[0] = R+S;
}

void _blowfly_dmeasure (double *lik, double *y, double *x, double *p, int give_log,
			int *obsindex, int *stateindex, int *parindex, int *covindex,
			int ncovars, double *covars, double t) {
  double size = 1.0/SIGMAY/SIGMAY;
  double prob = size/(size+N[0]);
  *lik = dnbinom(Y,size,prob,give_log);
}

void _blowfly_rmeasure (double *y, double *x, double *p, 
			  int *obsindex, int *stateindex, int *parindex, int *covindex,
			  int ncovars, double *covars, double t) {
  double size = 1.0/SIGMAY/SIGMAY;
  double prob = size/(size+N[0]);
  Y = rnbinom(size,prob);
}

#undef N
#undef R
#undef S
#undef E
#undef EPS

#undef P
#undef NZERO
#undef DELTA
#undef SIGMAP
#undef SIGMAD
#undef SIGMAY

void _blowfly_simulator_one (double *x, const double *p, 
			     const int *stateindex, const int *parindex, const int *covindex,
			     int covdim, const double *covar, 
			     double t, double dt) {
  _blowfly_simulator(x, p, stateindex, parindex, covindex, covdim, covar, t, dt, 14);
}

void _blowfly_simulator_two (double *x, const double *p, 
			     const int *stateindex, const int *parindex, const int *covindex,
			     int covdim, const double *covar, 
			     double t, double dt) {
  _blowfly_simulator(x, p, stateindex, parindex, covindex, covdim, covar, t, dt, 7);
}

