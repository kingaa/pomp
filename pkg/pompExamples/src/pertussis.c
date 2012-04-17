// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

static double expit (double x) {
  return 1.0/(1.0 + exp(-x));
}

static double logit (double x) {
  return log(x/(1-x));
}

static inline double rgammawn (double sigma, double dt) {
  double sigmasq;
  sigmasq = sigma*sigma;
  return (sigmasq > 0) ? rgamma(dt/sigmasq,sigmasq) : dt;
}

static double term_time (double t) {
  double day = 365.0 * (t - floor(t));
  int tt = (day >= 7.0 && day < 100.0) 
    || (day >= 116.0 && day < 200.0) 
    || (day >= 252.0 && day < 300.0) 
    || (day >= 308.0 && day < 356.0);
  return (2.0*((double) tt)-1.0);
}

#define BIRTHRATE        (p[parindex[0]]) // per capita birthrate
#define DEATHRATE        (p[parindex[1]]) // per capita deathrate
#define MEANBETA         (p[parindex[2]]) // mean contact rate
#define AMPLBETA         (p[parindex[3]]) // amplitude of contact rate
#define IMPORTS          (p[parindex[4]]) // annual influx of infectives
#define SIGMA            (p[parindex[5]]) // 1/LP
#define GAMMA            (p[parindex[6]]) // 1/IP
#define ALPHA1           (p[parindex[7]]) // 1/mean residence time in R1
#define ALPHA2           (p[parindex[8]]) // 1/mean residence time in V
#define ALPHA_RATIO      (p[parindex[9]]) // 1/mean residence time in R2
#define REPORT_PROB      (p[parindex[10]]) // reporting probability
#define BOOST_PROB       (p[parindex[11]]) // probability of immune boosting
#define POLAR_PHI        (p[parindex[12]]) // probability of failing to develop immunity
#define PVAC             (p[parindex[13]]) // vaccine coverage (probability)
#define FOI_MOD          (p[parindex[14]]) // force of infection modulator
#define NOISE_SIGMA      (p[parindex[15]]) // white noise intensity
#define POP              (p[parindex[16]]) // population size
#define TAU              (p[parindex[17]]) // reporting SD

#define SUS              (x[stateindex[0]]) // susceptibles
#define EXP              (x[stateindex[1]]) // exposed
#define INF              (x[stateindex[2]]) // infectious
#define R1               (x[stateindex[3]]) // 1st recovered class
#define R2               (x[stateindex[4]]) // 1st recovered class
#define VAC              (x[stateindex[5]]) // vaccinated
#define CASE             (x[stateindex[6]]) // cumulative cases
#define W                (x[stateindex[7]]) // cumulative noise
#define ERR              (x[stateindex[8]]) // errors
#define SIMPOP           (x[stateindex[9]]) // simulated population size

#define DSDT             (f[stateindex[0]]) // susceptibles
#define DEDT             (f[stateindex[1]]) // exposed
#define DIDT             (f[stateindex[2]]) // infectious
#define DR1DT            (f[stateindex[3]]) // 1st recovered class
#define DR2DT            (f[stateindex[4]]) // 1st recovered class
#define DVDT             (f[stateindex[5]]) // vaccinated
#define DCASEDT          (f[stateindex[6]]) // cumulative cases
#define DWDT             (f[stateindex[7]]) // cumulative noise
#define DERRDT           (f[stateindex[8]]) // errors
#define DSIMPOPDT        (f[stateindex[9]]) // simulated population size

#define REPORTS (y[obsindex[0]])

void pertussis_sveirr_EM (double *x, double *p, 
			  int *stateindex, int *parindex, int *covindex, 
			  int covdim, double *covar, double t, double dt)
{
  int ntrans = 17;
  double rate[ntrans];			     // transition rates
  //  int ntrans = sizeof(rate)/sizeof(rate[0]); // number of transitions
  double trans[ntrans];			     // transition numbers
  double beta, alpha3;
  double dW;		// white noise process
  double foi;		// force of infection
  void (*reulermultinom)(int, double, double *, double, double *);

  reulermultinom = (void (*)(int,double,double*,double,double*)) R_GetCCallable("pomp","reulermultinom");
  
  // untransform the parameters
  alpha3 = ALPHA1*ALPHA_RATIO;

  if (!(R_FINITE(BIRTHRATE))) return;
  if (!(R_FINITE(DEATHRATE))) return;
  if (!(R_FINITE(MEANBETA))) return;
  if (!(R_FINITE(AMPLBETA))) return;
  if (!(R_FINITE(IMPORTS))) return;
  if (!(R_FINITE(SIGMA))) return;
  if (!(R_FINITE(GAMMA))) return;
  if (!(R_FINITE(ALPHA1))) return;
  if (!(R_FINITE(ALPHA2))) return;
  if (!(R_FINITE(alpha3))) return;
  if (!(R_FINITE(BOOST_PROB))) return;
  if (!(R_FINITE(POLAR_PHI))) return;
  if (!(R_FINITE(PVAC))) return;
  if (!(R_FINITE(FOI_MOD))) return;
  if (!(R_FINITE(NOISE_SIGMA))) return;

  if (!(R_FINITE(SUS))) return;
  if (!(R_FINITE(EXP))) return;
  if (!(R_FINITE(INF))) return;
  if (!(R_FINITE(R1))) return;
  if (!(R_FINITE(R2))) return;
  if (!(R_FINITE(VAC))) return;
  if (!(R_FINITE(CASE))) return;
  if (!(R_FINITE(W))) return;

  beta = MEANBETA*(1+AMPLBETA*term_time(t));
  dW = rgammawn(NOISE_SIGMA,dt);
  if (!(R_FINITE(dW))) {
    ERR += 1.0e18;
    return;
  }

  foi = (IMPORTS+INF*beta*dW/dt)/POP; // force of infection

  // compute the transition rates
  // rates into S & V
  rate[0] = (1-PVAC)*BIRTHRATE*POP;   // birth into susceptible class
  rate[1] = PVAC*BIRTHRATE*POP;	      // birth into vaccinated class
  // rates out of S
  rate[2] = foi;		      // force of infection
  rate[3] = DEATHRATE;		      // death from susceptible class
  // rates out of E
  rate[4] = SIGMA;		      // termination of latent period
  rate[5] = DEATHRATE;		      // death from exposed class
  // rates out of I
  rate[6] = POLAR_PHI*GAMMA;	      // recovery to susceptible class
  rate[7] = (1-POLAR_PHI)*GAMMA;      // recovery to immune class
  rate[8] = DEATHRATE;		      // death from infectious class
  // rates out of R1
  rate[9] = ALPHA1;		      // waning of immunity
  rate[10] = DEATHRATE;		      // death from recovered class
  // rates out of R2
  rate[11] = BOOST_PROB*FOI_MOD*foi;	 // boosting of immunity
  rate[12] = (1-BOOST_PROB)*FOI_MOD*foi; // infection
  rate[13] = alpha3;			 // waning of R2 immunity
  rate[14] = DEATHRATE;			 // death from R2 
  // rates out of V
  rate[15] = ALPHA2;		      // waning of vaccine immunity
  rate[16] = DEATHRATE;		      // death from vaccinated class

  // compute the transition numbers
  trans[0] = rpois(rate[0]*dt);	// births are Poisson
  trans[1] = rpois(rate[1]*dt);
  reulermultinom(2,SUS, &rate[2],dt,&trans[2]);
  reulermultinom(2,EXP, &rate[4],dt,&trans[4]);
  reulermultinom(3,INF, &rate[6],dt,&trans[6]);
  reulermultinom(2,R1,  &rate[9],dt,&trans[9]);
  reulermultinom(4,R2, &rate[11],dt,&trans[11]);
  reulermultinom(2,VAC,&rate[15],dt,&trans[15]);

  // balance the equations
  SUS += trans[0] + trans[6] + trans[13] - trans[2] - trans[3];
  EXP += trans[2] + trans[12] - trans[4] - trans[5];
  INF += trans[4] - trans[6] - trans[7] - trans[8];
  R1 += trans[7] + trans[11] - trans[9] - trans[10];
  R2 += trans[9] + trans[15] - trans[11] - trans[12] - trans[13] - trans[14];
  VAC += trans[1] - trans[15] - trans[16];
  CASE += trans[4];
  W += (NOISE_SIGMA!=0) ? (dW-dt)/NOISE_SIGMA : 0.0; // mean zero, variance = dt

  SIMPOP = SUS + EXP + INF + R1 + R2 + VAC;
  // keep track of errors
  // if (SUS < 0.0) ERR += 1.0e0;
  // if (EXP < 0.0) ERR += 1.0e3;
  // if (INF < 0.0) ERR += 1.0e6;
  // if (R1 < 0.0)  ERR += 1.0e9;
  // if (R2 < 0.0)  ERR += 1.0e12;
  // if (VAC < 0.0) ERR += 1.0e15;
  ERR += 1.0;
}

void pertussis_sveirr_skel (double *f, double *x, double *p, 
			    int *stateindex, int *parindex, int *covindex, 
			    int covdim, double *covar, double t)
{
  int ntrans = 17;
  double rate[ntrans];			     // transition rates
  double alpha3;
  double beta;
  double foi;

  alpha3 = ALPHA1*ALPHA_RATIO;
  beta = MEANBETA*(1+AMPLBETA*term_time(t));
  foi = (IMPORTS+INF*beta)/POP; // force of infection

  // compute the transition rates
  // rates into S & V
  rate[0] = (1-PVAC)*BIRTHRATE*POP;   // birth into susceptible class
  rate[1] = PVAC*BIRTHRATE*POP;	      // birth into vaccinated class
  // rates out of S
  rate[2] = foi*SUS;		      // force of infection
  rate[3] = DEATHRATE*SUS;	      // death from susceptible class
  // rates out of E
  rate[4] = SIGMA*EXP;		      // termination of latent period
  rate[5] = DEATHRATE*EXP;	      // death from exposed class
  // rates out of I
  rate[6] = POLAR_PHI*GAMMA*INF; // recovery to susceptible class
  rate[7] = (1-POLAR_PHI)*GAMMA*INF;  // recovery to immune class
  rate[8] = DEATHRATE*INF;	      // death from infectious class
  // rates out of R1
  rate[9] = ALPHA1*R1;		      // waning of immunity
  rate[10] = DEATHRATE*R1;	      // death from recovered class
  // rates out of R2
  rate[11] = BOOST_PROB*FOI_MOD*foi*R2;	    // boosting of immunity
  rate[12] = (1-BOOST_PROB)*FOI_MOD*foi*R2; // infection
  rate[13] = alpha3*R2;			    // waning of R2 immunity
  rate[14] = DEATHRATE*R2;		    // death from R2 
  // rates out of V
  rate[15] = ALPHA2*VAC;	// waning of vaccine immunity
  rate[16] = DEATHRATE*VAC;	// death from vaccinated class

  // balance the equations
  DSDT = rate[0] + rate[6] + rate[13] - rate[2] - rate[3];
  DEDT = rate[2] + rate[12] - rate[4] - rate[5];
  DIDT = rate[4] - rate[6] - rate[7] - rate[8];
  DR1DT = rate[7] + rate[11] - rate[9] - rate[10];
  DR2DT = rate[9] + rate[15] - rate[11] - rate[12] - rate[13] - rate[14];
  DVDT = rate[1] - rate[15] - rate[16];

  DCASEDT = rate[4];	// roughly the weekly number of new cases
  DWDT = 0;
  DERRDT = 0;
  DSIMPOPDT = DSDT + DEDT + DIDT + DR1DT + DR2DT + DVDT;

}

void rounded_normal_dmeasure (double *lik, double *y, double *x, double *p, int give_log,
			      int *obsindex, int *stateindex, int *parindex, int *covindex,
			      int ncovars, double *covars, double t)
{
  double tol = 1e-20;
  double mn = REPORT_PROB*CASE;
  double sd = TAU*CASE;

  if (R_FINITE(sd) && R_FINITE(mn)) {
    if (ISNA(y[0]))
      lik[0] = 1.0;
    else if (y[0] > 0.0) 
      lik[0] = pnorm(y[0]+0.5,mn,sd+tol,1,0)-pnorm(y[0]-0.5,mn,sd+tol,1,0)+tol;
    else
      lik[0] = pnorm(y[0]+0.5,mn,sd+tol,1,0)+tol;
  } else {
    lik[0] = tol;
  }
  if (give_log) lik[0] = log(lik[0]);
}

void rounded_normal_rmeasure (double *y, double *x, double *p, 
			      int *obsindex, int *stateindex, int *parindex, int *covindex,
			      int ncovars, double *covars, double t)
{
  double tol = 1e-20;
  double mn = REPORT_PROB*CASE;
  double sd = TAU*CASE;

  if (R_FINITE(sd) && R_FINITE(mn)) {
    y[0] = rnorm(mn,sd+tol);
    if (y[0] > 0.0) 
      y[0] = nearbyint(y[0]);
    else
      y[0] = 0.0;
  } else {
    y[0] = R_NaReal;
  }
}

void negbin_dmeasure (double *lik, double *y, double *x, double *p, int give_log,
		      int *obsindex, int *stateindex, int *parindex, int *covindex,
		      int ncovars, double *covars, double t)
{
  double mu = REPORT_PROB*CASE;
  double size = 1.0/TAU;
  double prob = 1.0/(1.0+mu/size);
  double tol = 1e-20;

  if (ISNA(y[0])) {
    lik[0] = 1.0;
  } else if ((!(R_FINITE(prob))) || (!(R_FINITE(size)))) {
    *lik = tol;
  } else {
    *lik = dnbinom(y[0],size,prob,0)+tol;
  }
  if (give_log) *lik = log(*lik);
}

void negbin_rmeasure (double *y, double *x, double *p, 
		      int *obsindex, int *stateindex, int *parindex, int *covindex,
		      int ncovars, double *covars, double t)
{
  double mu = REPORT_PROB*CASE;
  double size = 1.0/TAU;
  double prob = 1.0/(1.0+mu/size);
  
  if ((!(R_FINITE(prob))) || (!(R_FINITE(size)))) {
    y[0] = R_NaReal;
  } else {
    y[0] = rnbinom(size,prob);
  }
}
