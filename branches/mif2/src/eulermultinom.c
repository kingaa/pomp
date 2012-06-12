#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

#include "pomp_internal.h"

void reulermultinom (int m, double size, double *rate, double dt, double *trans) {
  double p = 0.0;
  int j, k;
  if ((size < 0.0) || (dt < 0.0) || (floor(size+0.5) != size)) {
    for (k = 0; k < m; k++) trans[k] = R_NaN;
    return;
  }  
  for (k = 0; k < m; k++) {
    if (rate[k] < 0.0) {
      for (j = 0; j < m; j++) trans[j] = R_NaN;
      return;
    }
    p += rate[k]; // total event rate
  }
  if (p > 0.0) {
    size = rbinom(size,1-exp(-p*dt)); // total number of events
    if (!(R_FINITE(size)))
      warning("reulermultinom: result of binomial draw is not finite");
    m -= 1;
    for (k = 0; k < m; k++) {
      if (rate[k] > p) p = rate[k];
      trans[k] = ((size > 0) && (p > 0)) ? rbinom(size,rate[k]/p) : 0;
      if (!(R_FINITE(size)&&R_FINITE(p)&&R_FINITE(rate[k])&&R_FINITE(trans[k])))
	warning("reulermultinom: result of binomial draw is not finite");
      size -= trans[k];
      p -= rate[k];
    }
    trans[m] = size;
  } else {
    for (k = 0; k < m; k++) trans[k] = 0.0;
  }
}

// probability density of Euler-multinomial transitions
double deulermultinom (int m, double size, double *rate, double dt, double *trans, int give_log) {
  double p = 0.0;
  double n = 0.0;
  double ff = 0.0;
  int k;
  if ((dt < 0.0) || (size < 0.0) || (floor(size+0.5) != size)) {
    warning("NaNs produced");
    return R_NaN;
  }
  for (k = 0; k < m; k++) {
    if (rate[k] < 0.0) {
      warning("NaNs produced");
      return R_NaN;
    }
    if (trans[k] < 0.0) {
      ff = (give_log) ? R_NegInf: 0.0;
      return ff;
    }
    p += rate[k]; // total event rate
    n += trans[k]; // total number of events
  }
  if (n > size) {
    ff = (give_log) ? R_NegInf: 0.0;
    return ff;
  }
  ff = dbinom(n,size,1-exp(-p*dt),1); // total number of events
  m -= 1;
  for (k = 0; k < m; k++) {
    if ((n > 0) && (p > 0)) {
      if (rate[k] > p) p = rate[k];
      ff += dbinom(trans[k],n,rate[k]/p,1);
    } else if (trans[k] > 0.0) {
      ff = R_NegInf;
      return ff;
    }
    n -= trans[k];
    p -= rate[k];
  }
  ff = (give_log) ? ff : exp(ff);
  return ff;
}

void reulermultinom_multi (int *n, int *ntrans, double *size, double *rate, double *dt, double *trans) {
  int k;
  int m = *ntrans;
  for (k = 0; k < *n; k++) {
    reulermultinom(*ntrans,*size,rate,*dt,trans);
    trans += m;
  }
}

void deulermultinom_multi (int *n, int *ntrans, double *size, double *rate, double *dt, double *trans, int *give_log, double *f) {
  int k;
  int m = *ntrans;
  for (k = 0; k < *n; k++) {
    f[k] = deulermultinom(*ntrans,*size,rate,*dt,trans,*give_log);
    trans += m;
  }
}

SEXP R_Euler_Multinom (SEXP n, SEXP size, SEXP rate, SEXP dt) {
  int nprotect = 0;
  int ntrans = length(rate);
  int dim[2];
  SEXP nn, x;
  if (length(size)>1)
    warning("reulermultinom: only the first element of 'size' is meaningful");
  if (length(dt)>1)
    warning("reulermultinom: only the first element of 'dt' is meaningful");
  PROTECT(nn = AS_INTEGER(n)); nprotect++;
  dim[0] = ntrans;
  dim[1] = INTEGER(nn)[0];
  PROTECT(x = makearray(2,dim)); nprotect++;
  setrownames(x,GET_NAMES(rate),2);
  GetRNGstate();
  reulermultinom_multi(INTEGER(nn),&ntrans,REAL(size),REAL(rate),REAL(dt),REAL(x));
  PutRNGstate();
  UNPROTECT(nprotect);
  return x;
}

SEXP D_Euler_Multinom (SEXP x, SEXP size, SEXP rate, SEXP dt, SEXP log) {
  int nprotect = 0;
  int ntrans = length(rate);
  int *dim, n;
  SEXP f;
  dim = INTEGER(GET_DIM(x));
  if (dim[0] != ntrans)
    error("deulermultinom: 'NROW(x)' should match length of 'rate'");
  n = dim[1];
  if (length(size)>1)
    warning("deulermultinom: only the first element of 'size' is meaningful");
  if (length(dt)>1)
    warning("deulermultinom: only the first element of 'dt' is meaningful");
  PROTECT(f = NEW_NUMERIC(n)); nprotect++;
  deulermultinom_multi(&n,&ntrans,REAL(size),REAL(rate),REAL(dt),REAL(x),INTEGER(log),REAL(f));
  UNPROTECT(nprotect);
  return f;
}

// This function draws a random increment of a gamma whitenoise process.
// This will have expectation=dt and variance=(sigma^2*dt)
// If dW = rgammawn(sigma,dt), then 
// mu dW/dt is a candidate for a random rate process within an
// Euler-multinomial context, i.e., 
// E[mu*dW] = mu*dt and Var[mu*dW] = mu*sigma^2*dt
SEXP R_GammaWN (SEXP n, SEXP sigma, SEXP deltat) {
  int nprotect = 0;
  int k, nval, nsig, ndt;
  double *x, *sig, *dt;
  SEXP ans;
  PROTECT(n = AS_INTEGER(n)); nprotect++;
  nval = INTEGER(n)[0];
  PROTECT(sigma = AS_NUMERIC(sigma)); nprotect++;
  nsig = LENGTH(sigma);
  sig = REAL(sigma);
  PROTECT(deltat = AS_NUMERIC(deltat)); nprotect++;
  ndt = LENGTH(deltat);
  dt = REAL(deltat);
  PROTECT(ans = NEW_NUMERIC(nval)); nprotect++;
  x = REAL(ans);
  for (k = 0; k < nval; k++) 
    x[k] = rgammawn(sig[k%nsig],dt[k%ndt]);
  UNPROTECT(nprotect);
  return ans;
}

