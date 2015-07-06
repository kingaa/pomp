// dear emacs, please treat this as -*- C++ -*-

#include <Rmath.h>

#include "pomp.h"

#define B       (p[parindex[0]]) // egg-laying rate
#define MUL     (p[parindex[1]]) // larva survival
#define MUA     (p[parindex[2]]) // adult survival
#define CEA     (p[parindex[3]]) // adult-on-egg cannibalism
#define CEL     (p[parindex[4]]) // larva-on-egg cannibalism
#define CPA     (p[parindex[5]]) // adult-on-pupa cannibalism

#define L           (x[stateindex[0]]) // larvae
#define P           (x[stateindex[1]]) // pupae
#define A           (x[stateindex[2]]) // adults

#define LPRIME      (f[stateindex[0]]) // larvae
#define PPRIME      (f[stateindex[1]]) // pupae
#define APRIME      (f[stateindex[2]]) // adults

#define LOBS        (y[obsindex[0]]) // observed larvae 
#define POBS        (y[obsindex[1]]) // observed pupae
#define AOBS        (y[obsindex[2]]) // observed adults

void _lpa_original_skeleton (double *f, double *x, const double *p, 
			     const int *stateindex, const int *parindex, const int *covindex,
			     int covdim, const double *covar, double t) 
{
  double repro = B*exp(-CEL*L-CEA*A);
  double lsurv = 1.0-MUL;
  double psurv = exp(-CPA*A);
  double asurv = 1.0-MUA;
  
  LPRIME = A*repro;
  PPRIME = L*lsurv;
  APRIME = P*psurv+A*asurv;
  
}

// Poisson-binomial LPA model (Dennis et al. 2001)
void _lpa_poisson_binomial_simulator (double *x, const double *p, 
				      const int *stateindex, const int *parindex, const int *covindex,
				      int covdim, const double *covar, 
				      double t, double dt)
{
  double larvae, pupae, adults;

  double repro = B*exp(-CEL*L-CEA*A);
  double lsurv = 1.0-MUL;
  double psurv = exp(-CPA*A);
  double asurv = 1.0-MUA;
  
  larvae = rpois(A*repro);
  pupae = rbinom(L,lsurv);
  adults = rbinom(P,psurv)+rbinom(A,asurv);

  L = larvae;
  P = pupae;
  A = adults;
}

void _lpa_errorless_dmeasure (double *lik, double *y, double *x, double *p, int give_log,
			      int *obsindex, int *stateindex, int *parindex, int *covindex,
			      int ncovars, double *covars, double t) {
  double like;
  if ((fabs(L-LOBS)<0.5)&&
      (fabs(P-POBS)<0.5)&&
      (fabs(A-AOBS)<0.5))
    like = 1.0;
  else
    like = 0.0;

  *lik = (give_log) ? log(like) : like;
}

void _lpa_errorless_rmeasure (double *y, double *x, double *p, 
			      int *obsindex, int *stateindex, int *parindex, int *covindex,
			      int ncovars, double *covars, double t) {
  LOBS = nearbyint(L);
  POBS = nearbyint(P);
  AOBS = nearbyint(A);
}

void _lpa_adults_only_poisson_dmeasure (double *lik, double *y, double *x, double *p, int give_log,
					int *obsindex, int *stateindex, int *parindex, int *covindex,
					int ncovars, double *covars, double t) {
  *lik = dpois(AOBS,A,give_log);
}

void _lpa_adults_only_poisson_rmeasure (double *y, double *x, double *p, 
					int *obsindex, int *stateindex, int *parindex, int *covindex,
					int ncovars, double *covars, double t) {
  LOBS = 0;
  POBS = 0;
  AOBS = rpois(A);
}
