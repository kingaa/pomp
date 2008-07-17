#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

#include "pomp_internal.h"

// simulate Euler-multinomial transitions
void reulermultinom (int m, double size, double *rate, double dt, double *trans) {
  double p = 0.0;
  int j, k;
  if ((size < 0.0) || (dt < 0.0) || (floor(size+0.5) != size)) {
    for (k = 0; k < m; k++)
      trans[k] = R_NaN;
    return;
  }  
  for (k = 0; k < m; k++) {
    if (rate[k] < 0.0) {
      for (j = 0; j < m; j++)
	trans[j] = R_NaN;
      return;
    }
    p += rate[k]; // total event rate
  }
  GetRNGstate();
  size = rbinom(size,1-exp(-p*dt)); // total number of events
  if (!(R_FINITE(size)))
    Rprintf("reulermultinom infelicity 1: %lg %lg %lg\n",size,p,dt);
  m -= 1;
  for (k = 0; k < m; k++) {
    trans[k] = ((size > 0) && (p > 0)) ? rbinom(size,rate[k]/p) : 0;
    if (!(R_FINITE(size)&&R_FINITE(p)&&R_FINITE(rate[k])&&R_FINITE(trans[k])))
      Rprintf("reulermultinom infelicity %d: %lg %lg %lg %lg\n",k,size,rate[k]/p,trans[k],dt);
    size -= trans[k];
    p -= rate[k];
  }
  trans[m] = size;
  PutRNGstate();
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
  for (k = 0; k < *n; k++)
    reulermultinom(*ntrans,*size,rate,*dt,&trans[k*m]);
}

void deulermultinom_multi (int *n, int *ntrans, double *size, double *rate, double *dt, double *trans, int *give_log, double *f) {
  int k;
  int m = *ntrans;
  for (k = 0; k < *n; k++)
    f[k] = deulermultinom(*ntrans,*size,rate,*dt,&trans[k*m],*give_log);
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
  PROTECT(x =  makearray(2,dim)); nprotect++;
  setrownames(x,GET_NAMES(rate),2);
  reulermultinom_multi(INTEGER(nn),&ntrans,REAL(size),REAL(rate),REAL(dt),REAL(x));
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
  PROTECT(f =  NEW_NUMERIC(n)); nprotect++;
  deulermultinom_multi(&n,&ntrans,REAL(size),REAL(rate),REAL(dt),REAL(x),INTEGER(log),REAL(f));
  UNPROTECT(nprotect);
  return f;
}
