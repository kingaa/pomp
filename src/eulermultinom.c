#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

// simulate Euler-multinomial transitions
void reulermultinom (int m, double size, double *rate, double dt, double *trans) {
  double p = 0.0;
  int j, k;
  if ((size < 0.0) || (dt < 0.0)) {
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
void deulermultinom (int m, double size, double *rate, double *trans, double dt, int give_log, double *f) {
  double p = 0.0;
  double n = 0.0;
  double ff = 0.0;
  int k;
  if ((dt < 0.0) || (size < 0.0)) {
    warning("NaNs produced");
    *f = R_NaN;
    return;
  }
  for (k = 0; k < m; k++) {
    if (rate[k] < 0.0) {
      warning("NaNs produced");
      *f = R_NaN;
      return;
    }
    if (trans[k] < 0.0) {
      *f = (give_log) ? R_NegInf: 0.0;
      return;
    }
    p += rate[k]; // total event rate
    n += trans[k]; // total number of events
  }
  if (n > size) {
    *f = (give_log) ? R_NegInf: 0.0;
    return;
  }
  ff = dbinom(n,size,1-exp(-p*dt),1); // total number of events
  m -= 1;
  for (k = 0; k < m; k++) {
    if ((n > 0) && (p > 0)) {
      ff += dbinom(trans[k],n,rate[k]/p,1);
    } else if (trans[k] > 0.0) {
      ff = R_NegInf;
      return;
    }
    n -= trans[k];
    p -= rate[k];
  }
  *f = (give_log) ? ff : exp(ff);
}

void reulermultinom_multi (int *n, int *ntrans, double *size, double *rate, double *dt, double *trans) {
  int k;
  int m = *ntrans;
  for (k = 0; k < *n; k++)
    reulermultinom(*ntrans,*size,rate,*dt,&trans[k*m]);
}

void deulermultinom_multi (int *n, int *ntrans, double *size, double *rate, double *trans, 
			   double *dt, int *give_log, double *f) {
  int k;
  int m = *ntrans;
  for (k = 0; k < *n; k++)
    deulermultinom(*ntrans,*size,rate,&trans[k*m],*dt,*give_log,&f[k]);
}

