// -*- mode: C++ -*-

#include "pomp_internal.h"
#include <stdio.h>

static void pomp_acf (double *acf, double *x, int n, int nvars, int maxlag, int correlation);

SEXP probe_acf (SEXP x, SEXP lag_max, SEXP corr) {
  int nprotect = 0;
  SEXP acf, acf_names, cacf;
  SEXP X_names, Dim;
  int maxlag, correlation, nvars, n;
  int j, k, l;
  double *p, *p1;
  char tmp[BUFSIZ];

  maxlag = *(INTEGER(AS_INTEGER(lag_max)));    // maximum lag
  correlation = *(INTEGER(AS_INTEGER(corr))); // correlation, or covariance?

  PROTECT(Dim = GET_DIM(x)); nprotect++;
  nvars = INTEGER(Dim)[0]; 	// nvars = # of variables
  n = INTEGER(Dim)[1];		// n = # of observations
  // depending on how this function is called, it may be necessary to duplicate x in the next line
  PROTECT(x = AS_NUMERIC(x)); nprotect++; 
  PROTECT(X_names = GET_ROWNAMES(GET_DIMNAMES(x))); nprotect++;
   
  PROTECT(acf = NEW_NUMERIC((maxlag+1)*nvars)); nprotect++;

  pomp_acf(REAL(acf),REAL(x),n,nvars,maxlag,correlation);
  
  if (correlation) {
    PROTECT(cacf = NEW_NUMERIC(maxlag*nvars)); nprotect++;
    PROTECT(acf_names = NEW_STRING(maxlag*nvars)); nprotect++;
    for (j = 0, l = 0, p = REAL(cacf), p1 = REAL(acf)+1; j < nvars; j++, p1++) {
      for (k = 1; k <= maxlag; k++, p++, p1++, l++) {
	*p = *p1;		// copy correlations from acf into cacf
	snprintf(tmp,BUFSIZ,"acf.%ld.%s",k,CHARACTER_DATA(STRING_ELT(X_names,j)));
	SET_STRING_ELT(acf_names,l,mkChar(tmp));
      }
    }
    SET_NAMES(cacf,acf_names);
  } else {
    PROTECT(acf_names = NEW_STRING((maxlag+1)*nvars)); nprotect++;
    for (j = 0, l = 0; j < nvars; j++) {
      for (k = 0; k <= maxlag; k++, l++) {
	snprintf(tmp,BUFSIZ,"acf.%ld.%s",k,CHARACTER_DATA(STRING_ELT(X_names,j)));
	SET_STRING_ELT(acf_names,l,mkChar(tmp));
      }
    }
    SET_NAMES(acf,acf_names);
  }

  UNPROTECT(nprotect);
  if (correlation) return(cacf); 
  else return(acf);
}

// vectorized routine for ACF calculation
// thanks to Simon N. Wood for the original version of this code
// modifications due to AAK
static void pomp_acf (double *acf, double *x, int n, int nvars, int maxlag, int correlation) {
  double sum, *p, *p0, *p1, *p2;
  int j, k, lag, ct;

  // first center each row
  for (j = 0, p = x; j < nvars; j++, p++) {
    for (k = 0, p0 = p, sum = 0, ct = 0; k < n; p0 += nvars, k++) {
      if (R_FINITE(*p0)) {
	sum += *p0;
	ct++;
      }
    }
    if (ct < 1) error("series %ld has no data",j);
    sum /= ((double) ct);	// mean of x[j,]
    for (k = 0, p0 = p; k < n; p0 += nvars, k++)
      if (R_FINITE(*p0)) *p0 -= sum;
  }

  // compute covariances
  for (j = 0, p0 = x, p = acf; j < nvars; j++, p0++) { // loop over series
    for (lag = 0; lag <= maxlag; lag++, p++) { // loop over lags
      for (k = 0, ct = 0, sum = 0, p1 = p0, p2 = p0+lag*nvars; k < n-lag; k++, p1 += nvars, p2 += nvars)
  	if (R_FINITE(*p1) && R_FINITE(*p2)) {
  	  sum += (*p1)*(*p2);
  	  ct++;
  	}
      *p = (ct > 0) ? sum/((double) ct) : R_NaReal;
      ////      *p = (ct > 0) ? sum/((double) n) : R_NaReal; // strangely, this is apparently what R's 'acf' function does
    }
  }
  
  // convert to correlations if desired
  if (correlation) 
    for (j = 0, p = acf; j < nvars; j++, p += maxlag+1)
      for (lag = 0, p0 = p, sum = *p; lag <= maxlag; lag++, p0++) 
	*p0 /= sum;

}
