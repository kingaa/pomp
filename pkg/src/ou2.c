// -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

#include "pomp_internal.h"

// prototypes

void _ou2_normal_rmeasure (double *y, double *x, double *p, 
			   int *obsindex, int *stateindex, int *parindex, int *covindex,
			   int ncovar, double *covar, double t);
void _ou2_normal_dmeasure (double *lik, double *y, double *x, double *p, int give_log, 
			   int *obsindex, int *stateindex, int *parindex, int *covindex,
			   int covdim, double *covar, double t);
void _ou2_adv (double *x, double *xstart, double *par, double *times, int *n, int *parindex);
void _ou2_pdf (double *d, double *X, double *par, double *times, int *n, int *parindex, int *give_log);
static void sim_ou2 (double *x,
		     double alpha1, double alpha2, double alpha3, double alpha4, 
		     double sigma1, double sigma2, double sigma3);
static double dens_ou2 (double *x1, double *x2,
			double alpha1, double alpha2, double alpha3, double alpha4, 
			double sigma1, double sigma2, double sigma3, int give_log);

#define ALPHA1     (pp[parindex[0]])
#define ALPHA2     (pp[parindex[1]])
#define ALPHA3     (pp[parindex[2]])
#define ALPHA4     (pp[parindex[3]])
#define SIGMA1     (pp[parindex[4]])
#define SIGMA2     (pp[parindex[5]])
#define SIGMA3     (pp[parindex[6]])

// advance the matrix of particles from times[0] to the other times given
// it is assumed that the times are consecutive (FIX THIS!)
void _ou2_adv (double *x, double *xstart, double *par, double *times, int *n, int *parindex)
{
  int nvar = n[0], npar = n[1], nrep = n[2], ntimes = n[3];
  double *xp, *pp;
  int i, j, k;
  GetRNGstate();       // initialize R's pseudorandom number generator
  for (j = 0; j < nrep; j++) {
    xp = &x[nvar*j];		// get address of j-th state vector
    for (i = 0; i < nvar; i++) xp[i] = xstart[i+nvar*j]; // copy xstart into the first slice of x
  }
  for (k = 1; k < ntimes; k++) {
    R_CheckUserInterrupt();
    for (j = 0; j < nrep; j++) {
      xp = &x[nvar*(j+nrep*k)];
      pp = &par[npar*j];
      for (i = 0; i < nvar; i++) xp[i] = x[i+nvar*(j+nrep*(k-1))];
      sim_ou2(xp,ALPHA1,ALPHA2,ALPHA3,ALPHA4,SIGMA1,SIGMA2,SIGMA3); // advance particle
    }
  }
  PutRNGstate();	  // finished with R's random number generator
}

// pdf of a single 2D OU transition
void _ou2_pdf (double *d, double *X, double *par, double *times, int *n, int *parindex, int *give_log)
{
  int nvar = n[0], npar = n[1], nrep = n[2], ntimes = n[3];
  double *x1, *x2, *pp;
  int j, k;
  for (k = 0; k < nrep; k++) {
    pp = &par[npar*k];	       // get address of k-th parameter vector
    x1 = &X[nvar*k];	     // get address of (0,0)-th state vector
    for (j = 1; j < ntimes; j++) {
      R_CheckUserInterrupt();
      x2 = &X[nvar*(k+nrep*j)]; // get address of (k,j)-th state vector
      d[k+nrep*(j-1)] = dens_ou2(x1,x2,ALPHA1,ALPHA2,ALPHA3,ALPHA4,SIGMA1,SIGMA2,SIGMA3,*give_log);
      x1 = x2;
    }
  }
}

#undef ALPHA1
#undef ALPHA2
#undef ALPHA3
#undef ALPHA4
#undef SIGMA1
#undef SIGMA2
#undef SIGMA3

#define X1    (x[stateindex[0]])
#define X2    (x[stateindex[1]])
#define TAU   (p[parindex[7]])
#define Y1    (y[obsindex[0]])
#define Y2    (y[obsindex[1]])

// bivariate normal measurement error density
void _ou2_normal_dmeasure (double *lik, double *y, double *x, double *p, int give_log, 
			   int *obsindex, int *stateindex, int *parindex, int *covindex,
			   int covdim, double *covar, double t) 
{
  double sd = fabs(TAU);
  double f = 0.0;
  f += (ISNA(Y1)) ? 0.0 : dnorm(Y1,X1,sd,1);
  f += (ISNA(Y2)) ? 0.0 : dnorm(Y2,X2,sd,1);
  *lik = (give_log) ? f : exp(f);
}

// bivariate normal measurement error simulator
void _ou2_normal_rmeasure (double *y, double *x, double *p, 
			   int *obsindex, int *stateindex, int *parindex, int *covindex,
			   int ncovar, double *covar, 
			   double t) 
{
  double sd = fabs(TAU);
  Y1 = rnorm(X1,sd);
  Y2 = rnorm(X2,sd);
}

#undef X1
#undef X2
#undef TAU
#undef Y1
#undef Y2

// simple 2D Ornstein-Uhlenbeck process simulation
static void sim_ou2 (double *x,
		     double alpha1, double alpha2, double alpha3, double alpha4, 
		     double sigma1, double sigma2, double sigma3)
{
  double eps[2], xnew[2];

  if (!(R_FINITE(x[0]))) return;
  if (!(R_FINITE(x[1]))) return;
  if (!(R_FINITE(alpha1))) return;
  if (!(R_FINITE(alpha2))) return;
  if (!(R_FINITE(alpha3))) return;
  if (!(R_FINITE(alpha4))) return;
  if (!(R_FINITE(sigma1))) return;
  if (!(R_FINITE(sigma2))) return;
  if (!(R_FINITE(sigma3))) return;

  eps[0] = rnorm(0,1);
  eps[1] = rnorm(0,1);

  xnew[0] = alpha1*x[0]+alpha3*x[1]+sigma1*eps[0];
  xnew[1] = alpha2*x[0]+alpha4*x[1]+sigma2*eps[0]+sigma3*eps[1];

  x[0] = xnew[0];
  x[1] = xnew[1];
}

// simple 2D Ornstein-Uhlenbeck process density
static double dens_ou2 (double *x1, double *x2,
			double alpha1, double alpha2, double alpha3, double alpha4, 
			double sigma1, double sigma2, double sigma3, int give_log)
{
  double eps[2], val;

  if (!(R_FINITE(x1[0]))) return R_NaReal;
  if (!(R_FINITE(x1[1]))) return R_NaReal;
  if (!(R_FINITE(x2[0]))) return R_NaReal;
  if (!(R_FINITE(x2[1]))) return R_NaReal;
  if (!(R_FINITE(alpha1))) return R_NaReal;
  if (!(R_FINITE(alpha2))) return R_NaReal;
  if (!(R_FINITE(alpha3))) return R_NaReal;
  if (!(R_FINITE(alpha4))) return R_NaReal;
  if (!(R_FINITE(sigma1))) return R_NaReal;
  if (!(R_FINITE(sigma2))) return R_NaReal;
  if (!(R_FINITE(sigma3))) return R_NaReal;

  // compute residuals
  eps[0] = x2[0]-alpha1*x1[0]-alpha3*x1[1];
  eps[1] = x2[1]-alpha2*x1[0]-alpha4*x1[1];

  // backsolve
  eps[0] /= sigma1;
  eps[1] -= sigma2*eps[0];
  eps[1] /= sigma3;

  val = dnorm(eps[0],0.0,1.0,1)+dnorm(eps[1],0.0,1.0,1)-log(sigma1)-log(sigma3);
  return ((give_log) ? val : exp(val));
}
