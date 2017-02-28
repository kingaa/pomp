// -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include "pomp_internal.h"

static void sim_ou2 (double *x1, double *x2,
		     double alpha1, double alpha2, double alpha3, double alpha4, 
		     double sigma1, double sigma2, double sigma3);
static double dens_ou2 (double x1, double x2, double z1, double z2,
			double alpha1, double alpha2, double alpha3, double alpha4, 
			double sigma1, double sigma2, double sigma3, int give_log);

// #define ALPHA1     (pp[0])
// #define ALPHA2     (pp[1])
// #define ALPHA3     (pp[2])
// #define ALPHA4     (pp[3])
// #define SIGMA1     (pp[4])
// #define SIGMA2     (pp[5])
// #define SIGMA3     (pp[6])

// advance the matrix of particles from times[0] to the other times given
// note that you cannot assume that you can only assume that times[k]-times[k-1]>=0
// i.e., you cannot assume that successive times are consecutive, nor can you assume that
// they are distinct.
// void ou2_adv (double *x, double *xstart, double *par, double *times, int *n)
// {
//   int nvar = n[0], npar = n[1], nrep = n[2], ntimes = n[3], incr;
//   double tnow, tgoal, dt = 1.0;
//   double *xp0, *xp1, *pp;
//   int i, j, k;

//   incr = nrep*nvar;

//   GetRNGstate();       // initialize R's pseudorandom number generator

//   for (j = 0; j < nrep; j++) {

//     R_CheckUserInterrupt();	// check for an interrupt signal

//     xp0 = &xstart[nvar*j];     // pointer to j-th starting state
//     xp1 = &x[nvar*j];	       // pointer to j-th state vector
//     pp = &par[npar*j];	       // pointer to j-th parameter vector

//     for (i = 0; i < nvar; i++) xp1[i] = xp0[i]; // copy xstart into the first slice of x

//     tnow = times[0];		// initial time

//     for (k = 1; k < ntimes; k++) { // loop over times
    
//       xp0 = xp1;
//       xp1 += incr;

//       for (i = 0; i < nvar; i++) xp1[i] = xp0[i]; // copy state vector

//       tgoal = times[k];

//       while (tnow < tgoal) {
// 	sim_ou2(&xp1[0],&xp1[1],ALPHA1,ALPHA2,ALPHA3,ALPHA4,SIGMA1,SIGMA2,SIGMA3); // advance state
// 	tnow += dt;		// advance time
//       }

//     }

//   }

//   PutRNGstate();	  // finished with R's random number generator
// }

// #undef ALPHA1
// #undef ALPHA2
// #undef ALPHA3
// #undef ALPHA4
// #undef SIGMA1
// #undef SIGMA2
// #undef SIGMA3

#define ALPHA1     (p[parindex[0]])
#define ALPHA2     (p[parindex[1]])
#define ALPHA3     (p[parindex[2]])
#define ALPHA4     (p[parindex[3]])
#define SIGMA1     (p[parindex[4]])
#define SIGMA2     (p[parindex[5]])
#define SIGMA3     (p[parindex[6]])
#define TAU        (p[parindex[7]])

#define X1    (stateindex[0])
#define X2    (stateindex[1])

#define Y1    (y[obsindex[0]])
#define Y2    (y[obsindex[1]])

// onestep simulator for use in 'discrete.time.sim' plug-in
void ou2_step (double *x, const double *p,
	       const int *stateindex, const int *parindex, const int *covindex,
	       int ncovars, const double *covars,
	       double t, double dt) 
{
  sim_ou2(&x[X1],&x[X2],ALPHA1,ALPHA2,ALPHA3,ALPHA4,SIGMA1,SIGMA2,SIGMA3);
}

// onestep transition probability density for use in 'onestep.dens' plug-in
// transition from x to z as time goes from t1 to t2
void ou2_pdf (double *f, 
	      double *x, double *z, double t1, double t2, const double *p, 
	      const int *stateindex, const int *parindex, const int *covindex,
	      int ncovars, const double *covars)
{
  if (t2-t1 != 1)
    errorcall(R_NilValue,"ou2_pdf error: transitions must be consecutive");
  f[0] = dens_ou2(x[X1],x[X2],z[X1],z[X2],ALPHA1,ALPHA2,ALPHA3,ALPHA4,SIGMA1,SIGMA2,SIGMA3,1);
}

void ou2_skel (double *f, double *x, double *p, 
	       int *stateindex, int *parindex, int *covindex, 
	       int ncovars, double *covars, double t)
{
  f[X1] = ALPHA1*x[X1]+ALPHA3*x[X2];
  f[X2] = ALPHA2*x[X1]+ALPHA4*x[X2];
}

// bivariate normal measurement error density
void ou2_dmeasure (double *lik, double *y, double *x, double *p, int give_log, 
		   int *obsindex, int *stateindex, int *parindex, int *covindex,
		   int covdim, double *covar, double t) 
{
  double sd = fabs(TAU);
  double f = 0.0;
  f += (ISNA(Y1)) ? 0.0 : dnorm(Y1,x[X1],sd,1);
  f += (ISNA(Y2)) ? 0.0 : dnorm(Y2,x[X2],sd,1);
  *lik = (give_log) ? f : exp(f);
}

// bivariate normal measurement error simulator
void ou2_rmeasure (double *y, double *x, double *p, 
		   int *obsindex, int *stateindex, int *parindex, int *covindex,
		   int ncovar, double *covar, 
		   double t) 
{
  double sd = fabs(TAU);
  Y1 = rnorm(x[X1],sd);
  Y2 = rnorm(x[X2],sd);
}

#undef ALPHA1
#undef ALPHA2
#undef ALPHA3
#undef ALPHA4
#undef SIGMA1
#undef SIGMA2
#undef SIGMA3
#undef TAU

#undef X1
#undef X2
#undef Y1
#undef Y2

// simple 2D Ornstein-Uhlenbeck process simulation
static void sim_ou2 (double *x1, double *x2,
		     double alpha1, double alpha2, double alpha3, double alpha4, 
		     double sigma1, double sigma2, double sigma3)
{
  double eps[2], xnew[2];

  eps[0] = rnorm(0,1);
  eps[1] = rnorm(0,1);

  xnew[0] = alpha1*(*x1)+alpha3*(*x2)+sigma1*eps[0];
  xnew[1] = alpha2*(*x1)+alpha4*(*x2)+sigma2*eps[0]+sigma3*eps[1];

  *x1 = xnew[0];
  *x2 = xnew[1];
}

// simple 2D Ornstein-Uhlenbeck process transition density
// transition (x1,x2) -> (z1,z2) in 1 unit of time
static double dens_ou2 (double x1, double x2, double z1, double z2,
			double alpha1, double alpha2, double alpha3, double alpha4, 
			double sigma1, double sigma2, double sigma3, int give_log)
{
  double eps[2], val;

  // compute residuals
  eps[0] = z1-alpha1*x1-alpha3*x2;
  eps[1] = z2-alpha2*x1-alpha4*x2;

  // backsolve
  eps[0] /= sigma1;
  eps[1] -= sigma2*eps[0];
  eps[1] /= sigma3;

  val = dnorm(eps[0],0.0,1.0,1)+dnorm(eps[1],0.0,1.0,1)-log(sigma1)-log(sigma3);
  return ((give_log) ? val : exp(val));
}
