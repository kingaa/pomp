#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

#include "pomp_internal.h"

static void reulermultinom_multi (int *n, int *ntrans, double *size, double *rate, double *dt, double *trans) {
  int k;
  int m = *ntrans;
  for (k = 0; k < *n; k++) {
    reulermultinom(*ntrans,*size,rate,*dt,trans);
    trans += m;
  }
}

static void deulermultinom_multi (int *n, int *ntrans, double *size, double *rate, double *dt, double *trans, int *give_log, double *f) {
  int k;
  int m = *ntrans;
  for (k = 0; k < *n; k++) {
    f[k] = deulermultinom(*ntrans,*size,rate,*dt,trans,*give_log);
    trans += m;
  }
}

SEXP R_Euler_Multinom (SEXP n, SEXP size, SEXP rate, SEXP dt) {
  int ntrans = length(rate);
  int dim[2];
  SEXP nn, x, nm;
  if (length(size)>1)
    warn("in 'reulermultinom': only the first element of 'size' is meaningful");
  if (length(dt)>1)
    warn("in 'reulermultinom': only the first element of 'dt' is meaningful");
  PROTECT(nn = AS_INTEGER(n));
  dim[0] = ntrans;
  dim[1] = INTEGER(nn)[0];
  PROTECT(x = makearray(2,dim));
  PROTECT(nm = GET_NAMES(rate));
  setrownames(x,nm,2);
  GetRNGstate();
  reulermultinom_multi(INTEGER(nn),&ntrans,REAL(size),REAL(rate),REAL(dt),REAL(x));
  PutRNGstate();
  UNPROTECT(3);
  return x;
}

SEXP D_Euler_Multinom (SEXP x, SEXP size, SEXP rate, SEXP dt, SEXP log) {
  int ntrans = length(rate);
  int *dim, n;
  SEXP f;
  dim = INTEGER(GET_DIM(x));
  if (dim[0] != ntrans)
    err("NROW('x') should match length of 'rate'");
  n = dim[1];
  if (length(size)>1)
    warn("in 'deulermultinom': only the first element of 'size' is meaningful");
  if (length(dt)>1)
    warn("in 'deulermultinom': only the first element of 'dt' is meaningful");
  PROTECT(f = NEW_NUMERIC(n));
  deulermultinom_multi(&n,&ntrans,REAL(size),REAL(rate),REAL(dt),REAL(x),INTEGER(log),REAL(f));
  UNPROTECT(1);
  return f;
}

// This function draws a random increment of a gamma whitenoise process.
// This will have expectation=dt and variance=(sigma^2*dt)
// If dW = rgammawn(sigma,dt), then 
// mu dW/dt is a candidate for a random rate process within an
// Euler-multinomial context, i.e., 
// E[mu*dW] = mu*dt and Var[mu*dW] = mu*sigma^2*dt
SEXP R_GammaWN (SEXP n, SEXP sigma, SEXP deltat) {
  int k, nval, nsig, ndt;
  double *x, *sig, *dt;
  SEXP ans;
  PROTECT(n = AS_INTEGER(n));
  nval = INTEGER(n)[0];
  PROTECT(sigma = AS_NUMERIC(sigma));
  nsig = LENGTH(sigma);
  sig = REAL(sigma);
  PROTECT(deltat = AS_NUMERIC(deltat));
  ndt = LENGTH(deltat);
  dt = REAL(deltat);
  PROTECT(ans = NEW_NUMERIC(nval));
  x = REAL(ans);
  GetRNGstate();
  for (k = 0; k < nval; k++) {
    x[k] = rgammawn(sig[k%nsig],dt[k%ndt]);
  }
  PutRNGstate();
  UNPROTECT(4);
  return ans;
}

// This function draws a random variate from a beta-binomial distribution.
// The binomial has size n and probability P ~ Beta(mu,theta).
// That is, E[P] = mu and Var[P] = 
SEXP R_BetaBinom (SEXP n, SEXP size, SEXP prob, SEXP theta) {
  int k, nval, ns, np, nt;
  double *X, *S, *P, *T;
  SEXP ans;
  PROTECT(n = AS_INTEGER(n)); nval = INTEGER(n)[0];
  PROTECT(size = AS_NUMERIC(size)); ns = LENGTH(size); S = REAL(size);
  PROTECT(prob = AS_NUMERIC(prob)); np = LENGTH(prob); P = REAL(prob);
  PROTECT(theta = AS_NUMERIC(theta)); nt = LENGTH(theta); T = REAL(theta);
  PROTECT(ans = NEW_NUMERIC(nval)); X = REAL(ans);
  GetRNGstate();
  for (k = 0; k < nval; k++) {
    X[k] = rbetabinom(S[k%ns],P[k%np],T[k%nt]);
  }
  PutRNGstate();
  UNPROTECT(5);
  return ans;
}

SEXP D_BetaBinom (SEXP x, SEXP size, SEXP prob, SEXP theta, SEXP log) {
  int k, n, nx, ns, np, nt;
  double *F, *X, *S, *P, *T;
  SEXP f;
  PROTECT(x = AS_NUMERIC(x)); nx = LENGTH(x); X = REAL(x);
  PROTECT(size = AS_NUMERIC(size)); ns = LENGTH(size); S = REAL(size);
  PROTECT(prob = AS_NUMERIC(prob)); np = LENGTH(prob); P = REAL(prob);
  PROTECT(theta = AS_NUMERIC(theta)); nt = LENGTH(theta); T = REAL(theta);
  PROTECT(log = AS_INTEGER(log));
  n = (nx > ns) ? nx : ns;
  n = (n > np) ? n : np;
  n = (n > nt) ? n : nt;
  PROTECT(f = NEW_NUMERIC(n)); F = REAL(f);
  for (k = 0; k < n; k++) {
    F[k] = dbetabinom(X[k%nx],S[k%ns],P[k%np],T[k%nt],INTEGER(log)[0]);
  }
  UNPROTECT(6);
  return f;
}
