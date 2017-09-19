// -*- mode: C++ -*-

#include "pomp_internal.h"
#include <stdio.h>

// vectorized routine for ACF calculation
// thanks to Simon N. Wood for the original version of this code
// modifications due to AAK
// note that the behavior of this ACF is slightly different from that of R's 'acf' function
// we first center the series and then compute means of products
static void pomp_acf_compute (double *acf, double *x, int n, int nvars, int *lags, int nlag) {
  double xx, *p, *p0, *p1, *p2;
  int i, j, k, lag, ct;

  // first center each row
  for (j = 0, p = x; j < nvars; j++, p++) {
    for (k = 0, p0 = p, xx = 0, ct = 0; k < n; p0 += nvars, k++) {
      if (R_FINITE(*p0)) {
	xx += *p0;
	ct++;
      }
    }
    if (ct < 1) errorcall(R_NilValue,"series %ld has no data",j+1);
    xx /= ct;			// mean of x[j,]
    for (k = 0, p0 = p; k < n; p0 += nvars, k++)
      if (R_FINITE(*p0)) *p0 -= xx;
  }

  // compute covariances
  for (j = 0, p0 = x, p = acf; j < nvars; j++, p0++) { // loop over series
    for (i = 0; i < nlag; i++, p++) { // loop over lags
      lag = lags[i];		      // i-th lag
      for (k = 0, ct = 0, xx = 0, p1 = p0, p2 = p0+lag*nvars; k < n-lag; k++, p1 += nvars, p2 += nvars)
  	if (R_FINITE(*p1) && R_FINITE(*p2)) {
  	  xx += (*p1)*(*p2);
  	  ct++;
  	}
      *p = (ct > 0) ? xx/ct : R_NaReal;
    }
  }
  
}

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
  if (ct < 1) errorcall(R_NilValue,"series 1 has no data");
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
  if (ct < 1) errorcall(R_NilValue,"series 2 has no data");
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

SEXP probe_acf (SEXP x, SEXP lags, SEXP corr) {
  int nprotect = 0;
  SEXP ans, ans_names;
  SEXP X, X_names;
  int nlag, correlation, nvars, n;
  int j, k, l;
  double *p, *p1, *cov;
  int *lag;
  char tmp[BUFSIZ], *nm;

  nlag = LENGTH(lags);			      // number of lags
  PROTECT(lags = AS_INTEGER(lags)); nprotect++;
  lag = INTEGER(lags);
  correlation = *(INTEGER(AS_INTEGER(corr))); // correlation, or covariance?

  nvars = INTEGER(GET_DIM(x))[0]; 	// nvars = # of variables
  n = INTEGER(GET_DIM(x))[1];		// n = # of observations

  PROTECT(X = duplicate(AS_NUMERIC(x))); nprotect++; 
  PROTECT(X_names = GET_ROWNAMES(GET_DIMNAMES(x))); nprotect++;
   
  PROTECT(ans = NEW_NUMERIC(nlag*nvars)); nprotect++;

  pomp_acf_compute(REAL(ans),REAL(X),n,nvars,lag,nlag);

  if (correlation) {
    l = 0;
    cov = (double *) R_alloc(nvars,sizeof(double));
    pomp_acf_compute(cov,REAL(X),n,nvars,&l,1); // compute lag-0 covariance
    for (j = 0, p = REAL(ans), p1 = cov; j < nvars; j++, p1++)
      for (k = 0; k < nlag; k++, p++)
	*p /= *p1;
  }
  
  PROTECT(ans_names = NEW_STRING(nlag*nvars)); nprotect++;
  for (j = 0, l = 0; j < nvars; j++) {
    for (k = 0; k < nlag; k++, l++) {
      nm = (char *) CHAR(STRING_ELT(X_names,j));
      snprintf(tmp,BUFSIZ,"acf.%d.%s",lag[k],nm);
      SET_STRING_ELT(ans_names,l,mkChar(tmp));
    }
  }
  SET_NAMES(ans,ans_names);

  UNPROTECT(nprotect);
  return ans;
}

SEXP probe_ccf (SEXP x, SEXP y, SEXP lags, SEXP corr) {
  int nprotect = 0;
  SEXP ans, ans_names;
  SEXP X, Y;
  double cov[2], xx, *p;
  int nlag, n, correlation;
  int k;
  char tmp[BUFSIZ];
  
  nlag = LENGTH(lags);
  PROTECT(lags = AS_INTEGER(lags)); nprotect++;
  correlation = *(INTEGER(AS_INTEGER(corr))); // correlation, or covariance?

  n = LENGTH(x);		// n = # of observations
  if (n != LENGTH(y))
    errorcall(R_NilValue,"'x' and 'y' must have equal lengths");

  PROTECT(X = duplicate(AS_NUMERIC(x))); nprotect++; 
  PROTECT(Y = duplicate(AS_NUMERIC(y))); nprotect++; 
   
  PROTECT(ans = NEW_NUMERIC(nlag)); nprotect++;

  pomp_ccf_compute(REAL(ans),REAL(X),REAL(Y),n,INTEGER(lags),nlag);
  
  if (correlation) {
    k = 0;
    pomp_acf_compute(&cov[0],REAL(X),n,1,&k,1); // compute lag-0 covariance of x
    pomp_acf_compute(&cov[1],REAL(Y),n,1,&k,1); // compute lag-0 covariance of y
    xx = sqrt(cov[0]*cov[1]);
    for (k = 0, p = REAL(ans); k < nlag; k++, p++) *p /= xx; // convert to correlation
  }
  
  PROTECT(ans_names = NEW_STRING(nlag)); nprotect++;
  for (k = 0; k < nlag; k++) {
    snprintf(tmp,BUFSIZ,"ccf.%d",INTEGER(lags)[k]);
    SET_STRING_ELT(ans_names,k,mkChar(tmp));
  }
  SET_NAMES(ans,ans_names);

  UNPROTECT(nprotect);
  return ans;
}

