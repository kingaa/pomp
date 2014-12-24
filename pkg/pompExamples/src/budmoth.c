// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <pomp.h>

#define ALPHA      (p[parindex[0]])
#define SIGALPHA   (p[parindex[1]])
#define GAM        (p[parindex[2]])
#define LAMBDA     (p[parindex[3]])
#define SIGLAMBDA  (p[parindex[4]])
#define GEE        (p[parindex[5]])
#define DELTA      (p[parindex[6]])
#define AEY        (p[parindex[7]])
#define SIGAEY     (p[parindex[8]])
#define DUBYA      (p[parindex[9]])
#define BETA0      (p[parindex[10]])
#define BETA1      (p[parindex[11]])
#define YEW        (p[parindex[12]])
#define SIGQOBS    (p[parindex[13]])
#define SIGNOBS    (p[parindex[14]])
#define SIGSOBS    (p[parindex[15]])

#define ALPHASTATE     (x[stateindex[0]])
#define LAMBDASTATE    (x[stateindex[1]])
#define ASTATE         (x[stateindex[2]])
#define QSTATE         (x[stateindex[3]])
#define NSTATE         (x[stateindex[4]])
#define SSTATE         (x[stateindex[5]])

#define QOBS   (y[obsindex[0]])
#define NOBS   (y[obsindex[1]])
#define SOBS   (y[obsindex[2]])

void budmoth_map (double *x, double *p, 
		  int *stateindex, int *parindex, int *covindex, 
		  int covdim, double *covar, double t, double dt)
{
  double Q, N, S;
  double tol = 1e-6;
  double sig2lambda;

  // check the discrete-time assumption
  if (dt != 1.0) 
    error("'delta.t'=",dt,"!=1 violates an assumption of this model");

  sig2lambda = SIGLAMBDA*SIGLAMBDA/LAMBDA;

  ALPHASTATE = expit(rnorm(logit(ALPHA),SIGALPHA));
  LAMBDASTATE = (SIGLAMBDA != 0) ? rgamma(LAMBDA/sig2lambda,sig2lambda) : LAMBDA;
  ASTATE = rlnorm(log(AEY),SIGAEY);

  Q = (1-ALPHASTATE)*(GAM/(GAM+NSTATE))+ALPHASTATE*QSTATE;
  N = LAMBDASTATE*NSTATE*(1-SSTATE)*exp(-GEE*NSTATE-DELTA*(1-QSTATE));
  S = 1-tol-exp(-(ASTATE*SSTATE*NSTATE+2*tol)/(1+DUBYA*ASTATE*SSTATE*NSTATE));

  QSTATE = Q;
  NSTATE = N;
  SSTATE = S;

}

void budmoth_skeleton (double *f, double *x, double *p, 
		       int *stateindex, int *parindex, int *covindex, 
		       int ncovars, double *covars, double t) 
{
  double Q, N, S;
  double tol = 1e-6;

  f[stateindex[0]] = ALPHA;	// ALPHA equation
  f[stateindex[1]] = LAMBDA;	// LAMBDA equation
  f[stateindex[2]] = AEY;	// A equation
  f[stateindex[3]] = (1-ALPHA)*(GAM/(GAM+NSTATE))+ALPHA*QSTATE; // Q equation
  f[stateindex[4]] = LAMBDA*NSTATE*(1-SSTATE)*exp(-GEE*NSTATE-DELTA*(1-QSTATE)); // N equation
  f[stateindex[5]] = 1-tol-exp(-(ASTATE*SSTATE*NSTATE+2*tol)/(1+DUBYA*ASTATE*SSTATE*NSTATE)); // S equation

}

void budmoth_rmeasure (double *y, double *x, double *p, 
		       int *obsindex, int *stateindex, int *parindex, int *covindex,
		       int ncovars, double *covars, double t) {
  double meanlogQ = log(BETA0+BETA1*QSTATE);
  double meanlogN = log(NSTATE);
  double meanlogitS = logit(YEW*SSTATE);

  if (!(R_FINITE(meanlogQ)) || !(R_FINITE(SIGQOBS))) {
    QOBS = R_NaN;
  } else {
    QOBS = rlnorm(meanlogQ,SIGQOBS);
  }

  if (!(R_FINITE(meanlogN = log(NSTATE))) || !(R_FINITE(SIGNOBS))) {
    NOBS = R_NaN;
  } else {
    NOBS = rlnorm(meanlogN,SIGNOBS);
  }

  if (!(R_FINITE(meanlogitS)) || !(R_FINITE(SIGSOBS))) {
    SOBS = R_NaN;
  } else {
    SOBS = expit(rnorm(meanlogitS,SIGSOBS));
  }
  
}

void budmoth_dmeasure (double *lik, double *y, double *x, double *p, int give_log,
		       int *obsindex, int *stateindex, int *parindex, int *covindex,
		       int ncovars, double *covars, double t) {
  
  double meanlogQ = log(BETA0+BETA1*QSTATE);
  double meanlogN = log(NSTATE);
  double meanlogitS = logit(YEW*SSTATE);
  double f1, f2, f3;

  f1 = dlnorm(QOBS,meanlogQ,SIGQOBS,1);
  f2 = dlnorm(NOBS,meanlogN,SIGNOBS,1);
  f3 = dnorm(logit(SOBS),meanlogitS,SIGSOBS,1)-log(SOBS*(1-SOBS));
  
  *lik = (give_log) ? f1+f2+f3 : exp(f1+f2+f3);
}

#undef QOBS
#undef NOBS
#undef SOBS

#undef ALPHASTATE
#undef LAMBDASTATE
#undef ASTATE
#undef QSTATE
#undef NSTATE
#undef SSTATE

#define ALPHASTATE     (stateindex[0])
#define LAMBDASTATE    (stateindex[1])
#define ASTATE         (stateindex[2])
#define QSTATE         (stateindex[3])
#define NSTATE         (stateindex[4])
#define SSTATE         (stateindex[5])

void budmoth_density (double *f, double *x1, double *x2, double t1, double t2, const double *p, 
		      const int *stateindex, const int *parindex, const int *covindex,
		      int covdim, const double *covar)
{
  double noise;
  double sig2lambda;
  double Q, N, S;
  double tol = 1e-6;
  double f1, f2, f3;

  // check the discrete-time assumption
  if (t2-t1 != 1.0) 
    error("t2-t1=",t2-t1,"!=1 violates an assumption of this model");

  sig2lambda = SIGLAMBDA*SIGLAMBDA/LAMBDA;

  f1 = dnorm(logit(x2[ALPHASTATE]),logit(ALPHA),SIGALPHA/(x2[ALPHASTATE]*(1-x2[ALPHASTATE])),1);
  f2 = dgamma(x2[LAMBDASTATE],LAMBDA/sig2lambda,sig2lambda,1);
  f3 = dlnorm(x2[ASTATE],log(AEY),SIGAEY,1);
  *f = f1+f2+f3;
}

#undef ALPHASTATE
#undef LAMBDASTATE
#undef ASTATE
#undef QSTATE
#undef NSTATE
#undef SSTATE

#undef ALPHA
#undef SIGALPHA
#undef GAM
#undef LAMBDA
#undef SIGLAMBDA
#undef GEE
#undef DELTA
#undef AEY
#undef SIGAEY
#undef DUBYA
#undef BETA0
#undef BETA1
#undef YEW
#undef SIGQOBS
#undef SIGNOBS
#undef SIGSOBS
