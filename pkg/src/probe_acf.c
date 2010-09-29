// -*- mode: C++ -*-

#include "pomp_internal.h"
#include <stdio.h>

static void pomp_acf (double *acf, double *x, int n, int nvars, int maxlag, int correlation);

SEXP probe_acf (SEXP x, SEXP lag_max, SEXP corr) {
  int nprotect = 0;
  SEXP acf, acf_names, cacf;
  SEXP X, X_names;
  int maxlag, correlation, nvars, n;
  int j, k, l;
  double *p, *p1;
  char tmp[BUFSIZ], *nm;

  maxlag = *(INTEGER(AS_INTEGER(lag_max)));    // maximum lag
  correlation = *(INTEGER(AS_INTEGER(corr))); // correlation, or covariance?

  nvars = INTEGER(GET_DIM(x))[0]; 	// nvars = # of variables
  n = INTEGER(GET_DIM(x))[1];		// n = # of observations

  PROTECT(X = duplicate(AS_NUMERIC(x))); nprotect++; 
  PROTECT(X_names = GET_ROWNAMES(GET_DIMNAMES(x))); nprotect++;
   
  PROTECT(acf = NEW_NUMERIC((maxlag+1)*nvars)); nprotect++;

  pomp_acf(REAL(acf),REAL(X),n,nvars,maxlag,correlation);
  
  if (correlation) {
    PROTECT(cacf = NEW_NUMERIC(maxlag*nvars)); nprotect++;
    PROTECT(acf_names = NEW_STRING(maxlag*nvars)); nprotect++;
    for (j = 0, l = 0, p = REAL(cacf), p1 = REAL(acf)+1; j < nvars; j++, p1++) {
      for (k = 1; k <= maxlag; k++, p++, p1++, l++) {
	*p = *p1;		// copy correlations from acf into cacf
	nm = (char *) CHARACTER_DATA(STRING_ELT(X_names,j));
	snprintf(tmp,BUFSIZ,"acf.%d.%s",k,nm);
	SET_STRING_ELT(acf_names,l,mkChar(tmp));
      }
    }
    SET_NAMES(cacf,acf_names);
  } else {
    PROTECT(acf_names = NEW_STRING((maxlag+1)*nvars)); nprotect++;
    for (j = 0, l = 0; j < nvars; j++) {
      for (k = 0; k <= maxlag; k++, l++) {
	nm = (char *) CHARACTER_DATA(STRING_ELT(X_names,j));
	snprintf(tmp,BUFSIZ,"acf.%d.%s",k,nm);
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
// note that the behavior of this ACF is slightly different from that of R's 'acf' function
// we first center the series and then compute means of products
static void pomp_acf (double *acf, double *x, int n, int nvars, int maxlag, int correlation) {
  double xx, *p, *p0, *p1, *p2;
  int j, k, lag, ct;

  // first center each row
  for (j = 0, p = x; j < nvars; j++, p++) {
    for (k = 0, p0 = p, xx = 0, ct = 0; k < n; p0 += nvars, k++) {
      if (R_FINITE(*p0)) {
	xx += *p0;
	ct++;
      }
    }
    if (ct < 1) error("series %ld has no data",j);
    xx /= ct;			// mean of x[j,]
    for (k = 0, p0 = p; k < n; p0 += nvars, k++)
      if (R_FINITE(*p0)) *p0 -= xx;
  }

  // compute covariances
  for (j = 0, p0 = x, p = acf; j < nvars; j++, p0++) { // loop over series
    for (lag = 0; lag <= maxlag; lag++, p++) { // loop over lags
      for (k = 0, ct = 0, xx = 0, p1 = p0, p2 = p0+lag*nvars; k < n-lag; k++, p1 += nvars, p2 += nvars)
  	if (R_FINITE(*p1) && R_FINITE(*p2)) {
  	  xx += (*p1)*(*p2);
  	  ct++;
  	}
      *p = (ct > 0) ? xx/ct : R_NaReal;
    }
  }
  
  // convert to correlations if desired
  if (correlation) 
    for (j = 0, p = acf; j < nvars; j++, p += maxlag+1)
      for (lag = 0, p0 = p, xx = *p; lag <= maxlag; lag++, p0++) 
	*p0 /= xx;

}
