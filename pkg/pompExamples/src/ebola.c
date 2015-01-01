// SEIR Ebola model

#include <R.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <pomp.h>

// State variables
#define S (x[stateindex[0]]) // Susceptible
#define E(J) (x[stateindex[1] + (J)]) // Exposed
#define I (x[stateindex[2]]) // Infected
#define R (x[stateindex[3]]) // Removed
#define N_EI (x[stateindex[4]]) // Number of transitions from E to I
#define N_IR (x[stateindex[5]]) // Number of transitions from I to R

// Variations
#define DS (f[stateindex[0]]) // Susceptible
#define DE(J) (f[stateindex[1] + (J)]) // Exposed
#define DI (f[stateindex[2]]) // Infected
#define DR (f[stateindex[3]]) // Removed
#define DN_EI (f[stateindex[4]]) // Number of transitions from E to I
#define DN_IR (f[stateindex[5]]) // Number of transitions from I to R

// Parameters on the natural scale (all rates are per day)
#define N (p[parindex[0]]) // Population size
#define R0 (p[parindex[1]]) // Basic reproduction number
#define alpha (p[parindex[2]]) // Inverse of latency period
#define gamma (p[parindex[3]]) // Inverse of duration of infection
#define rho (p[parindex[4]]) // Reporting probability
#define k (p[parindex[5]]) // Reporting overdispersion
#define cfr (p[parindex[6]]) // Case-fatality ratio
#define IC(J) (p[parindex[7] + (J)]) // Initial conditions

// Parameters on the transformed scale (all rates are per day)
#define TN (pt[parindex[0]]) // Population size
#define TR0 (pt[parindex[1]]) // Basic reproduction number
#define Talpha (pt[parindex[2]]) // Inverse of latency period
#define Tgamma (pt[parindex[3]]) // Inverse of duration of infection
#define Trho (pt[parindex[4]]) // Reporting probability
#define Tk (pt[parindex[5]]) // Reporting overdispersion
#define Tcfr (pt[parindex[6]]) // Case-fatality ratio
#define TIC(J) (pt[parindex[7] + (J)]) // Initial conditions

// Observations
#define cases (y[obsindex[0]]) // Number of reported cases
#define deaths (y[obsindex[1]]) // Number of reported deaths

// Transforms the parameters to the transformed scale
void _ebola_par_untrans (double *pt, double *p, int *parindex){
  TN = log(N);
  TR0 = log(R0);
  Talpha = log(alpha);
  Tgamma = log(gamma);
  Trho = logit(rho);
  Tk = log(k);
  Tcfr = logit(cfr);
  to_log_barycentric(&(TIC(0)),&(IC(0)),4);
}

// Transforms the parameters to the natural scale
void _ebola_par_trans (double *pt, double *p, int *parindex){
  TN = exp(N);
  TR0 = exp(R0);
  Talpha = exp(alpha);
  Tgamma = exp(gamma);
  Trho = expit(rho);
  Tk = exp(k);
  Tcfr = expit(cfr);
  from_log_barycentric(&(TIC(0)),&(IC(0)),4);
}

//  Observation model: hierarchical model for cases and deaths
// p(R_t, D_t| C_t) = p(R_t | C_t) * p(D_t | C_t, R_t)
// p(R_t | C_t): Negative binomial with mean rho * C_t and dispersion parameter 1 / k
// p(D_t | C_t, R_t): Binomial B(R_t, cfr)
void _ebola_dObs (double *lik, double *y, double *x, double *p, int give_log,
		  int *obsindex, int *stateindex, int *parindex, int *covindex,
		  int ncovars, double *covars, double t) {
  double f;
  if (k > 0.0)
    f = dnbinom_mu(nearbyint(cases),1.0/k,rho*N_EI,1);
  // f += dnbinom_mu(nearbyint(deaths), 1.0 / k, rho * cfr * N_IR, 1);
  else
    f = dpois(nearbyint(cases),rho*N_EI,1);
  *lik = (give_log) ? f : exp(f);
}

// For least-squares trajectory-matching:
void _ebola_dObsLS (double *lik, double *y, double *x, double *p, int give_log,
		    int *obsindex, int *stateindex, int *parindex, int *covindex,
		    int ncovars, double *covars, double t) {
  double f;
  f = dnorm(cases,rho*N_EI,k,1);
  *lik = (give_log) ? f : exp(f);
}

void _ebola_rObs (double *y, double *x, double *p,
		  int *obsindex, int *stateindex, int *parindex, int *covindex,
		  int ncovars, double *covars, double t)
{
  if (k > 0) {
    cases = rnbinom_mu(1.0/k,rho*N_EI);
    deaths = rnbinom_mu(1.0/k,rho*cfr*N_IR);
  } else {
    cases = rpois(rho*N_EI);
    deaths = rpois(rho*cfr*N_IR);
  }
}

// For least-squares trajectory-matching:
void _ebola_rObsLS (double *y, double *x, double *p,
		    int *obsindex, int *stateindex, int *parindex, int *covindex,
		    int ncovars, double *covars, double t)
{
  cases = rnorm(rho*N_EI,k);
  deaths = NA_REAL;
}


// Process model
void _ebola_rSim (double *x, const double *p,
		  const int *stateindex, const int *parindex, const int *covindex,
		  int covdim, const double *covars,
		  double t, double dt) {

  // Retrieve user data in the pomp object
  int *(*get_pomp_userdata_int)(const char *);
  get_pomp_userdata_int = (int *(*)(const char *)) R_GetCCallable("pomp","get_pomp_userdata_int");
  int nstageE = *(get_pomp_userdata_int("nstageE")); // Number of stages in the E class

  // Other parameters
  double lambda, beta;
  beta = R0 * gamma; // Transmission rate
  lambda = beta * I / N; // Force of infection
  int i;

  // Transitions

  // From class S
  double transS = rbinom(S, 1.0 - exp(- lambda * dt)); // No of infections

  // From class E
  double transE[nstageE]; // No of transitions between classes E
  for(i = 0; i < nstageE; i++){
    transE[i] = rbinom(E(i), 1.0 - exp(- nstageE * alpha * dt));
  }

  // From class I
  double transI = rbinom(I, 1.0 - exp(- gamma * dt)); // No of transitions I->R

  // Balance the equations
  S -= transS;
  E(0) += transS - transE[0];
  for(i=1; i < nstageE; i++) {
    E(i) += transE[i-1] - transE[i];
  }
  I += transE[nstageE - 1] - transI;
  R += transI;
  N_EI += transE[nstageE - 1]; // No of transitions from E to I
  N_IR += transI; // No of transitions from I to R
}

// Continuous-time deterministic skeleton
void _ebola_skel (double *f, double *x, double *p,
		  int *stateindex, int *parindex, int *covindex,
		  int ncovars, double *covars, double t) {

  // Retrieve user data in the pomp object
  int *(*get_pomp_userdata_int)(const char *);
  get_pomp_userdata_int = (int *(*)(const char *)) R_GetCCallable("pomp","get_pomp_userdata_int");
  int nstageE = *(get_pomp_userdata_int("nstageE")); // Number of stages in the E class

  // Other parameters
  double lambda, beta;
  beta = R0 * gamma; // Transmission rate
  lambda = beta * I / N; // Force of infection
  int i;

  // Balance the equations
  DS = - lambda * S;
  DE(0) = lambda * S - nstageE * alpha * E(0);
  for(i=1; i < nstageE; i++) DE(i) = nstageE * alpha * (E(i - 1) -  E(i));
  DI = nstageE * alpha * E(nstageE - 1) - gamma * I;
  DR = gamma * I;
  DN_EI = nstageE * alpha * E(nstageE - 1);
  DN_IR = gamma * I;
}
