// -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include "pomp.h"

double expit (double x) {return 1.0/(1.0+exp(-x));}

#define TAU        (p[parindex[0]])
#define GAMMA      (p[parindex[1]])
#define EPS        (p[parindex[2]])
#define DELTA      (p[parindex[3]])
#define DELTA_I    (p[parindex[4]])
#define OMEGA      (p[parindex[5]])
#define SD_BETA    (p[parindex[6]])
#define BETATREND  (p[parindex[7]])
#define LOGBETA    (p[parindex[8]])
#define ALPHA      (p[parindex[9]])
#define RHO        (p[parindex[10]])
#define CLIN       (p[parindex[11]])
#define NBASIS     (p[parindex[12]])
#define NRSTAGE    (p[parindex[13]])

#define SUSCEP     (x[stateindex[0]])
#define INFECT     (x[stateindex[1]])
#define RSHORT     (x[stateindex[2]])
#define RLONG      (x[stateindex[3]])
#define DEATHS     (x[stateindex[4]])
#define NOISE      (x[stateindex[5]])
#define COUNT      (x[stateindex[6]])

#define POP        (covar[covindex[0]])
#define DPOPDT     (covar[covindex[1]])
#define SEASBASIS  (covar[covindex[2]])
#define TREND      (covar[covindex[3]])

#define DATADEATHS   (y[obsindex[0]])

void _cholmodel_norm_rmeasure (double *y, double *x, double *p, 
			       int *obsindex, int *stateindex, int *parindex, int *covindex,
			       int ncovars, double *covars, double t)
{
  double v, tol = 1.0e-18;
  v = DEATHS*exp(TAU);
  if ((COUNT > 0) || (!(R_FINITE(v)))) {
    DATADEATHS = R_NaReal;
  } else {
    DATADEATHS = rnorm(DEATHS,v+tol);
  }
}

void _cholmodel_norm_dmeasure (double *lik, double *y, double *x, double *p, int give_log,
			       int *obsindex, int *stateindex, int *parindex, int *covindex,
			       int ncovars, double *covars, double t)
{
  double v, tol = 1.0e-18;
  v = DEATHS*exp(TAU);
  if ((COUNT>0.0) || (!(R_FINITE(v)))) {
    *lik = tol;
  } else {
    *lik = dnorm(DATADEATHS,DEATHS,v+tol,0)+tol;
  }
  if (give_log) *lik = log(*lik);
}

#undef DATADEATHS

// two-path SIRS cholera model using SDEs
// exponent (alpha) on I/n
// only "severe" infections are infectious
// truncation is not used
// instead, particles with negative states are killed
void _cholmodel_one (double *x, const double *p, 
		const int *stateindex, const int *parindex, const int *covindex,
		int covdim, const double *covar, 
		double t, double dt)
{			   // implementation of the SIRS cholera model
  int nrstage = (int) NRSTAGE;
  int nbasis  = (int) NBASIS;
  double births;
  double infections;
  double sdeaths;
  double ideaths;
  double rsdeaths;
  double rldeaths[nrstage];
  double disease;
  double wanings;
  double passages[nrstage+1];
  double effI;
  double beta;
  double omega;
  double gamma;
  double rho;
  double clin;
  double eps; 
  double delta; 
  double deltaI;
  double sd_beta;
  double alpha;
  double dw;
  const double *pt;
  int j;

  if (!(R_FINITE(SUSCEP))) return;
  if (!(R_FINITE(INFECT))) return;
  if (!(R_FINITE(RSHORT))) return;
  for (pt = &RLONG, j = 0; j < nrstage; j++) {
    if (!(R_FINITE(pt[j]))) return;
  }

  gamma = exp(GAMMA);
  eps = exp(EPS)*NRSTAGE;
  rho = exp(RHO);
  delta = exp(DELTA);
  deltaI = exp(DELTA_I);
  clin = expit(CLIN);
  sd_beta = exp(SD_BETA);
  alpha = exp(ALPHA);

  if (!(R_FINITE(gamma))) return;
  if (!(R_FINITE(eps))) return;
  if (!(R_FINITE(rho))) return;
  if (!(R_FINITE(delta))) return;
  if (!(R_FINITE(deltaI))) return;
  if (!(R_FINITE(clin))) return;
  if (!(R_FINITE(sd_beta))) return;
  if (!(R_FINITE(alpha))) return;
  if (!(R_FINITE(BETATREND))) return;

  beta = exp(dot_product(nbasis,&SEASBASIS,&LOGBETA)+BETATREND*TREND);
  omega = exp(dot_product(nbasis,&SEASBASIS,&OMEGA));

  if (!(R_FINITE(beta))) return;
  if (!(R_FINITE(omega))) return;

  dw = rnorm(0,sqrt(dt));	// white noise

  effI = pow(INFECT/POP,alpha);
  births = DPOPDT + delta*POP;	// births

  passages[0] = gamma*INFECT;	// recovery
  ideaths = delta*INFECT;	// natural I deaths
  disease = deltaI*INFECT;	// disease death
  rsdeaths = delta*RSHORT;	// natural Rs deaths
  wanings = rho*RSHORT;		// loss of immunity

  for (pt = &RLONG, j = 0; j < nrstage; j++) {
    rldeaths[j] = delta*pt[j];	// natural R deaths
    passages[j+1] = eps*pt[j];	// passage to the next immunity class
  }

  infections = (omega+(beta+sd_beta*dw/dt)*effI)*SUSCEP; // infection
  sdeaths = delta*SUSCEP;	// natural S deaths
  if (infections > 0.0) {
    if ((infections+sdeaths)*dt > SUSCEP) { // too many leaving S class
      COUNT += 1.0e10;
      return;
    }
  } else {
    if ((-clin*infections+disease+ideaths+passages[0])*dt > INFECT) { // too many leaving I class
      COUNT += 1.0e5;
      return;
    }
    if ((-(1-clin)*infections+rsdeaths+wanings)*dt > RSHORT) { // too many leaving Rs class
      COUNT += 1;
      return;
    }
  }

  SUSCEP += (births - infections - sdeaths + passages[nrstage] + wanings)*dt;
  INFECT += (clin*infections - disease - ideaths - passages[0])*dt;
  RSHORT += ((1-clin)*infections - rsdeaths - wanings)*dt;
  for (j = 0; j < nrstage; j++) 
    (&RLONG)[j] += (passages[j] - passages[j+1] - rldeaths[j])*dt;
  DEATHS += disease*dt;		// cumulative deaths due to disease
  NOISE += dw;

}

#undef GAMMA    
#undef EPS      
#undef DELTA    
#undef DELTA_I  
#undef OMEGA    
#undef SD_BETA  
#undef BETATREND
#undef LOGBETA  
#undef ALPHA    
#undef RHO      
#undef CLIN     
#undef NBASIS   
#undef NRSTAGE  

#undef SUSCEP   
#undef INFECT   
#undef RSHORT   
#undef RLONG    
#undef DEATHS   
#undef NOISE    
#undef COUNT    

#undef POP      
#undef DPOPDT   
#undef SEASBASIS
