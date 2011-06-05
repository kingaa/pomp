// dear emacs, please treat this as -*- C++ -*-

#include <Rmath.h>

#include "pomp.h"

#define LOG_BETA    (p[parindex[0]]) // transmission rates
#define LOG_M       (p[parindex[1]]) // Poisson import rate
#define LOG_ALPHA   (p[parindex[2]]) // mixing exponent
#define LOG_RHO     (p[parindex[3]]) // under-reporting
#define PERIOD      (p[parindex[4]]) // period of seasonality
#define DEGREE      (p[parindex[5]]) // degree of B-splines
#define NBASIS      (p[parindex[6]]) // number of B-spline basis functions
#define LOG_SIGMA   (p[parindex[7]]) // noise intensity
#define LOG_OD      (p[parindex[8]]) // measurement overdispersion parameter

#define BIRTHS      (covar[covindex[0]]) // numbers of births

#define SUS         (x[stateindex[0]]) // susceptibles
#define INF         (x[stateindex[1]]) // infectives
#define THETA       (x[stateindex[2]]) // imports

#define REPORTS     (y[obsindex[0]]) // reported cases

#define NEWSUS      (f[stateindex[0]]) // susceptibles
#define NEWINF      (f[stateindex[1]]) // infectives
#define NEWTHETA    (f[stateindex[2]]) // imports

void _tsir_dmeasure (double *lik, double *y, double *x, double *p, int give_log,
			  int *obsindex, int *stateindex, int *parindex, int *covindex,
			  int ncovars, double *covars, double t) 
{
  double mu, size, prob;
  size = exp(-LOG_OD);		// od = 1/size
  mu = exp(LOG_RHO)*INF;	// mean
  prob = size/(size+mu);
  //  *lik = dbinom(REPORTS,nearbyint(INF),exp(LOG_RHO),give_log);
  *lik = dnbinom(REPORTS,size,prob,give_log);
}

void _tsir_rmeasure (double *y, double *x, double *p, 
			  int *obsindex, int *stateindex, int *parindex, int *covindex,
			  int ncovars, double *covars, double t) {
  double mu, size, prob;
  size = exp(-LOG_OD);		// od = 1/size
  mu = exp(LOG_RHO)*INF;	// mean
  prob = size/(size+mu);
  //  REPORTS = rbinom(nearbyint(INF),exp(LOG_RHO));
  REPORTS = rnbinom(size,prob);
}

// Ricker model with log-normal process noise
void _tsir_simulator (double *x, const double *p, 
		      const int *stateindex, const int *parindex, const int *covindex,
		      int covdim, const double *covar, 
		      double t, double dt)
{
  double beta;
  double em;
  double alpha;
  double sigma, eps;
  double lambda;
  int nbasis;
  int degree;

  nbasis = (int) NBASIS;
  degree = (int) DEGREE;

  em = exp(LOG_M);
  alpha = exp(LOG_ALPHA);
  sigma = exp(LOG_SIGMA);
  
  {
    double seasonality[nbasis];
    periodic_bspline_basis_eval(t,PERIOD,degree,nbasis,&seasonality[0]);
    beta = exp(dot_product(nbasis,&seasonality[0],&LOG_BETA));
  }

  if (sigma!=0) 
    eps = exp(rnorm(0,sigma)); 
  else
    eps = 1;

  THETA = rpois(em);			  // imports
  lambda = beta*pow(INF+THETA,alpha)*eps; // force of infection
  INF = rbinom(nearbyint(SUS),1-exp(-lambda));
  SUS += BIRTHS-INF;		 // susceptible balance
}

void _tsir_skeleton (double *f, double *x, const double *p, 
		     const int *stateindex, const int *parindex, const int *covindex,
		     int covdim, const double *covar, double t) 
{
  double beta;
  double em;
  double alpha;
  double lambda;
  int nbasis;
  int degree;

  nbasis = (int) NBASIS;
  degree = (int) DEGREE;

  em = exp(LOG_M);
  alpha = exp(LOG_ALPHA);

  {
    double seasonality[nbasis];
    periodic_bspline_basis_eval(t,PERIOD,degree,nbasis,&seasonality[0]);
    beta = exp(dot_product(nbasis,&seasonality[0],&LOG_BETA));
  }

  lambda = beta*pow(INF+em,alpha); // force of infection
  NEWINF = SUS*(1-exp(-lambda));
  NEWSUS = SUS+BIRTHS-NEWINF;	// susceptible balance
  NEWTHETA = em;		// imports

}
