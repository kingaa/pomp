// -*- mode: C++ -*-

#include "pomp_internal.h"
#include <stdio.h>

// vectorized routine for CCF calculation
// we first center the series and then compute means of products
static void pomp_ccf_compute (double *ccf, double *x, double *y, int n, int *lags, int nlag) {
  double xx, *p, *p1, *p2;
  int j, k, lag, ct;

  // first center x
  for (k = 0, xx = 0, ct = 0, p = x; k < n; k++, p++) {
    if (R_FINITE(*p)) {
      xx += *p;
      ct++;
    }
  }
  if (ct < 1) error("series 'x' has no data");
  xx /= ct;			// mean of x[j]
  for (k = 0, p = x; k < n; k++, p++)
    if (R_FINITE(*p)) *p -= xx;

  // now center y
  for (k = 0, xx = 0, ct = 0, p = y; k < n; k++, p++) {
    if (R_FINITE(*p)) {
      xx += *p;
      ct++;
    }
  }
  if (ct < 1) error("series 'y' has no data");
  xx /= ct;			// mean of y[j]
  for (k = 0, p = y; k < n; k++, p++)
    if (R_FINITE(*p)) *p -= xx;
    
  // compute covariances
  for (j = 0, p = ccf; j < nlag; j++, p++) { // loop over lags
    lag = lags[j];
    if (lag < 0) {
      for (k = 0, xx = 0, ct = 0, p1 = x-lag, p2 = y; k < n+lag; k++, p1++, p2++)
	if (R_FINITE(*p1) && R_FINITE(*p1)) {
	  xx += (*p1)*(*p2);
	  ct++;
	}
      *p = (ct > 0) ? xx/ct : R_NaReal;
    } else {
      for (k = 0, xx = 0, ct = 0, p1 = x, p2 = y+lag; k < n-lag; k++, p1++, p2++)
	if (R_FINITE(*p1) && R_FINITE(*p1)) {
	  xx += (*p1)*(*p2);
	  ct++;
	}
      *p = (ct > 0) ? xx/ct : R_NaReal;
    }
  }
  
}

SEXP probe_ccf (SEXP x, SEXP y, SEXP lags) {
  int nprotect = 0;
  SEXP ccf, ccf_names;
  SEXP X, Y;
  int nlag, n;
  int k;
  char tmp[BUFSIZ], *nm;
  
  nlag = LENGTH(lags);
  PROTECT(lags = AS_INTEGER(lags)); nprotect++;

  n = LENGTH(x);		// n = # of observations
  if (n != LENGTH(y))
    error("'x' and 'y' must have equal lengths");

  PROTECT(X = duplicate(AS_NUMERIC(x))); nprotect++; 
  PROTECT(Y = duplicate(AS_NUMERIC(y))); nprotect++; 
   
  PROTECT(ccf = NEW_NUMERIC(nlag)); nprotect++;

  pomp_ccf_compute(REAL(ccf),REAL(X),REAL(Y),n,INTEGER(lags),nlag);
  
  PROTECT(ccf_names = NEW_STRING(nlag)); nprotect++;
  for (k = 0; k < nlag; k++) {
    snprintf(tmp,BUFSIZ,"ccf.%d",INTEGER(lags)[k]);
    SET_STRING_ELT(ccf_names,k,mkChar(tmp));
  }
  SET_NAMES(ccf,ccf_names);

  UNPROTECT(nprotect);
  return ccf;
}

