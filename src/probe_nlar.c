// -*- mode: C++; -*-

#include "internal.h"
#include "pomp_mat.h"
#include <stdio.h>

static void poly_nlar_fit(double *beta, double *y, int n, int nterms, int *lag, int *power, double *X);

SEXP probe_nlar (SEXP x, SEXP lags, SEXP powers) {
  SEXP y, beta, beta_names;
  int n, nterms;
  int k;
  double *mm;
  char tmp[BUFSIZ];

  n = LENGTH(x);		// n = # of observations
  nterms = LENGTH(lags);

  PROTECT(y = duplicate(AS_NUMERIC(x)));
  PROTECT(beta = NEW_NUMERIC(nterms));
  
  mm = (double *) R_alloc(n*nterms,sizeof(double)); // storage for the model matrix
  poly_nlar_fit(REAL(beta),REAL(y),n,nterms,INTEGER(lags),INTEGER(powers),mm);
  
  PROTECT(beta_names = NEW_STRING(nterms));
  for (k = 0; k < nterms; k++) {
    snprintf(tmp,BUFSIZ,"nlar.%d^%d",INTEGER(lags)[k],INTEGER(powers)[k]);
    SET_STRING_ELT(beta_names,k,mkChar(tmp));
  }
  SET_NAMES(beta,beta_names);

  UNPROTECT(3);
  return beta;
}

// Code for polynomial auto-regression
// The original version of the following code is due to Simon N. Wood.
// Modifications by AAK.
static void poly_nlar_fit (double *beta, double *y, int n, 
			   int nterms, int *lag, int *power, double *X) {
  // 'x' is an n vector of data.
  // 'nterms' gives the number of terms on the rhs of the autoregression.
  // 'lag[i]' gives the lag of the ith term on the rhs.
  // 'power[i]' gives the power to which the ith rhs term should be raised.
  // 'beta' contains the ar coefficients

  int maxlag, nx, ny, ok, ct;
  int i, j, k;
  double xx, *yp;
  double obs1;

  // find maxlag
  for (maxlag = 0, i = 0; i < nterms; i++) 
    maxlag = (lag[i] > maxlag) ? lag[i] : maxlag;
  ny = n - maxlag;		// maximum response vector length

  // compute the series mean
  for (j = 0, ct = 0, xx = 0.0; j < n; j++) 
    if (R_FINITE(y[j])) { 
      xx += y[j];
      ct++;
    }
  xx /= ct; // series mean

  // center the whole series
  // also check to see if there is any variability in the predictors
  ok = 0;
  obs1 = R_NaReal;
  for (j = 0; j < n; j++) {
    if (R_FINITE(y[j])) {
      if (!ok && (j < ny)) { 	// j < ny means x[j] is a predictor
	if (!R_FINITE(obs1)) obs1 = y[j]; 
	else if (y[j] != obs1) ok = 1;
      } 
      y[j] -= xx;		// subtracting series mean
    } 
  }
  
  if (!ok) {			// data had no variability
    
    for (i = 0; i < nterms; i++) beta[i] = 0.0;

  } else {			// data not all the same
      
    double *Xp;
    int finite[ny];
  
    // test for NA rows in model matrix and response vector
    for (nx = 0, yp = y+maxlag, j = 0; j < ny; j++) {
      finite[j] = (R_FINITE(yp[j])) ? 1 : 0; // finite response?
      for (i = 0; i < nterms; i++)
	finite[j] = (R_FINITE(yp[j-lag[i]])) ? finite[j] : 0; // finite in model matrix row j, column i
      if (finite[j]) nx++;
    }
    // nx is now the number of non-NA rows in the model matrix

    // build the model matrix, omitting NA rows
    for (Xp = X, i = 0; i < nterms; i++) { // work through the terms
      for (j = 0; j < ny; j++) {
	if (finite[j]) {
	  xx = yp[j-lag[i]];	// current predictor
	  *Xp = xx;
	  for (k = 1; k < power[i]; k++) *Xp *= xx; // raise to the appropriate power
	  Xp++;
	}
      }
    }
    // X is now the nx by nterms model matrix

    // drop the NA rows from the response data
    for (i = 0, j = 0; j < ny; j++)
      if (finite[j]) yp[i++] = yp[j]; // keep this row
    // response vector is now length nx

    {
      double tau[nterms];
      int pivot[nterms];

      // first QR decompose the model matrix
      pomp_qr(X,nx,nterms,pivot,tau);
      // then solve R b = Q'y for b
      pomp_qrqy(yp,X,nx,tau,nx,1,nterms,"left","transpose"); // y <- Q'y 
      pomp_backsolve(X,nx,nterms,yp,1,"Upper","No transpose","Non-unit");   // y <- R^{-1} y
      
      // unpivot and store coefficients in beta
      for (i = 0; i < nterms; i++) beta[pivot[i]] = yp[i];

    }

  }

}
