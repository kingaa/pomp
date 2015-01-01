// -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include "pomp.h"

#define TAU        (p[parindex[0]])
#define GAMMA      (p[parindex[1]])
#define EPS        (p[parindex[2]])
#define DELTA      (p[parindex[3]])
#define DELTA_I    (p[parindex[4]])
#define LOGOMEGA   (p[parindex[5]])
#define SD_BETA    (p[parindex[6]])
#define BETATREND  (p[parindex[7]])
#define LOGBETA    (p[parindex[8]])
#define ALPHA      (p[parindex[9]])
#define RHO        (p[parindex[10]])
#define CLIN       (p[parindex[11]])
#define NBASIS     (p[parindex[12]])
#define NRSTAGE    (p[parindex[13]])
#define S0         (p[parindex[14]])
#define I0         (p[parindex[15]])
#define RS0        (p[parindex[16]])
#define RL0        (p[parindex[17]])

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

void _cholmodel_untrans (double *pt, double *p, int *parindex) 
{
  int k, nrstage = (int) NRSTAGE;
  pt[parindex[0]] = log(TAU);
  pt[parindex[1]] = log(GAMMA);
  pt[parindex[2]] = log(EPS);
  pt[parindex[3]] = log(DELTA);
  pt[parindex[4]] = log(DELTA_I);
  pt[parindex[6]] = log(SD_BETA);
  pt[parindex[9]] = log(ALPHA);
  pt[parindex[10]] = log(RHO);
  pt[parindex[11]] = logit(CLIN);

  to_log_barycentric(&pt[parindex[14]],&S0,3+nrstage);
}
 
void _cholmodel_trans (double *pt, double *p, int *parindex) 
{
  int k, nrstage = (int) NRSTAGE;
  pt[parindex[0]] = exp(TAU);
  pt[parindex[1]] = exp(GAMMA);
  pt[parindex[2]] = exp(EPS);
  pt[parindex[3]] = exp(DELTA);
  pt[parindex[4]] = exp(DELTA_I);
  pt[parindex[6]] = exp(SD_BETA);
  pt[parindex[9]] = exp(ALPHA);
  pt[parindex[10]] = exp(RHO);
  pt[parindex[11]] = expit(CLIN);

  from_log_barycentric(&pt[parindex[14]],&S0,3+nrstage);
}

void _cholmodel_norm_rmeasure (double *y, double *x, double *p, 
			       int *obsindex, int *stateindex, 
			       int *parindex, int *covindex,
			       int ncovars, double *covars, double t)
{
  double tol = 1.0e-18;
  if ((COUNT > 0) || (!(R_FINITE(DEATHS)))) {
    DATADEATHS = R_NaReal;
  } else { 
    DATADEATHS = rnbinom_mu(1 / (TAU*TAU), DEATHS + tol);
  }
}

void _cholmodel_norm_dmeasure (double *lik, double *y, double *x, 
			       double *p, int give_log,
			       int *obsindex, int *stateindex, 
			       int *parindex, int *covindex,
			       int ncovars, double *covars, double t)
{
  double tol = 1.0e-18;
  if ((COUNT>0.0) || (!(R_FINITE(DEATHS)))) {
    *lik = tol;
  } else {
    *lik = dnbinom_mu(DATADEATHS+tol, 1 / (TAU*TAU), DEATHS + tol, 0) + tol;
  }
  if (give_log) *lik = log(*lik);
}

#undef DATADEATHS


inline double out_flow(double init, double rate_by_dt)
{

  return init * (1.0 - exp( - rate_by_dt ));

}

// two-path SIRS cholera model using SDEs
// exponent (alpha) on I/n
// only "severe" infections are infectious
// truncation is not used
// instead, particles with negative states are killed
void _cholmodel_two (double *x, const double *p, 
		     const int *stateindex, const int *parindex, 
		     const int *covindex,
		     int covdim, const double *covar, 
		     double t, double dt)
{			   // implementation of the SIRS cholera model
  int nrstage = (int) NRSTAGE;
  int nbasis  = (int) NBASIS;

  double *pt_r, *pt_w;
  int j;

  double k_eps = NRSTAGE * EPS;
 
  double gamma_delta_deltaI = GAMMA + DELTA + DELTA_I;

  double keps_div_keps_delta = k_eps / (k_eps + DELTA);

  /* 
   * Preliminaries
   */

  double beta = exp(dot_product(nbasis,&SEASBASIS,&LOGBETA)+BETATREND*TREND);

  double omega = exp(dot_product(nbasis,&SEASBASIS,&LOGOMEGA));

  double dw = rgammawn(SD_BETA, dt);  // gamma noise, mean=dt, variance=(beta_sd^2 dt)

  double lambda = omega + (beta * dw/dt) * INFECT/POP;  // Time-dependent force of infection

  /* 
   * Out-flows from compartments using aboundances at previous time step
   */

  double S_o = out_flow(SUSCEP, (lambda + DELTA) * dt);
  double I_o = out_flow(INFECT, gamma_delta_deltaI * dt);
  double RSH_o = out_flow(RSHORT, (RHO + DELTA) * dt);

  double R_o[nrstage];
  for (pt_r = &RLONG, pt_w = R_o, j = 0; 
       j < nrstage; 
       j++, pt_r++, pt_w++) {

    *pt_w = out_flow(*pt_r, (k_eps + DELTA) * dt); 

  }

  double new_infect = S_o * (lambda / (lambda + DELTA));

  /*
   * Cholera deaths = (I_o * DELTA_I) / gamma_delta_deltaI are not reborn in SUSCEPT.
   * To offset this downward bias on total population, we add the average number of (observed) deaths:
   *
   * average_chol_death = tot_cholc_deaths / (n_months * n_step)  <----->  29.522 = 354275.0 / (600 * 200) 
   */
  double chol_deaths_dt = 29.523;
  double births = DPOPDT * dt + chol_deaths_dt + S_o * (DELTA / (lambda + DELTA)) + I_o * DELTA / gamma_delta_deltaI + RSH_o * (DELTA / (DELTA + RHO));
  for (pt_r = R_o, j = 0; 
       j < nrstage; 
       j++, pt_r++) {

    births += *pt_r * (DELTA / (k_eps + DELTA)); 

  } 

  /* 
   * Updating all compartments
   */

  SUSCEP += - S_o + births + R_o[nrstage-1] * keps_div_keps_delta + RSH_o * (RHO / (RHO + DELTA));

  INFECT += - I_o + CLIN * new_infect;

  RSHORT += - RSH_o + (1 - CLIN) * new_infect;

  for (pt_w = &RLONG, pt_r = R_o, j = 0; 
       j < nrstage; 
       j++, pt_r++, pt_w++) {

    if(j == 0)
      {

	*pt_w += - *pt_r + I_o * (GAMMA / gamma_delta_deltaI); 

      } else {

      *pt_w += - *pt_r + *(pt_r-1) * keps_div_keps_delta;

    } 

  }

  DEATHS += I_o * (DELTA_I / gamma_delta_deltaI);	// cumulative deaths due to disease

  NOISE += (dw-dt) / SD_BETA;

}

#undef GAMMA    
#undef EPS      
#undef DELTA    
#undef DELTA_I  
#undef LOGOMEGA    
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
