#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

#include "pomp_internal.h"

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
    warning("in 'reulermultinom': only the first element of 'size' is meaningful");
  if (length(dt)>1)
    warning("in 'reulermultinom': only the first element of 'dt' is meaningful");
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
    error("NROW('x') should match length of 'rate'");
  n = dim[1];
  if (length(size)>1)
    warning("in 'deulermultinom': only the first element of 'size' is meaningful");
  if (length(dt)>1)
    warning("in 'deulermultinom': only the first element of 'dt' is meaningful");
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
