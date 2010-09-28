// -*- mode: C++; -*-

#include "pomp_internal.h"
#include <stdio.h>

static void pomp_nlar(double *beta, double *y, int n, int nterms, int *lag, int *power);

SEXP probe_nlar (SEXP x, SEXP lags, SEXP powers) {
  int nprotect = 0;
  SEXP y, beta, beta_names;
  int n, nterms;
  int k;
  char tmp[BUFSIZ];

  n = LENGTH(x);		// n = # of observations
  nterms = LENGTH(lags);

  PROTECT(y = duplicate(AS_NUMERIC(x))); nprotect++; 
  PROTECT(beta = NEW_NUMERIC(nterms)); nprotect++;

  pomp_nlar(REAL(beta),REAL(y),n,nterms,INTEGER(lags),INTEGER(powers));
  
  PROTECT(beta_names = NEW_STRING(nterms)); nprotect++;
  for (k = 0; k < nterms; k++) {
    snprintf(tmp,BUFSIZ,"nlar.%ld^%ld",INTEGER(lags)[k],INTEGER(powers)[k]);
    SET_STRING_ELT(beta_names,k,mkChar(tmp));
  }
  SET_NAMES(beta,beta_names);

  UNPROTECT(nprotect);
  return beta;
}

// Code for polynomial auto-regression
// The original version of the following code is due to Simon N. Wood.
// Modifications by AAK.
static void pomp_nlar(double *beta, double *y, int n, 
		      int nterms, int *lag, int *power) {
  // 'x' is an n vector of data.
  // 'nterms' gives the number of terms on the rhs of the autoregression.
  // 'lag[i]' gives the lag of the ith term on the rhs.
  // 'power[i]' gives the power to which the ith rhs term should be raised.
  // 'beta' contains the ar coefficients

  int maxlag, nx, ny, ok, ct;
  int one = 1;
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
      
    double *X, *Xp;
    int finite[ny];
  
    // test for NA rows in model matrix and response vector
    for (nx = 0, yp = y+maxlag, j = 0; j < ny; j++) {
      finite[j] = (R_FINITE(yp[j])) ? 1 : 0; // finite response?
      for (i = 0; i < nterms; i++) {
	finite[j] = (R_FINITE(yp[j-lag[i]])) ? finite[j] : 0; // finite in model matrix row j, column i
      }
      if (finite[j]) nx++;
    }
    // nx is now the number of non-NA rows in the model matrix

    // build the model matrix, omitting NA rows
    X = (double *) Calloc(nx*nterms,double);
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
    for (i = 0, j = 0; j < ny; j++) {
      if (finite[j]) {		// keep this row
	yp[i++] = yp[j];
      }
    }	
    // response vector is now length nx

    {
      double tau[nterms], b[nterms];
      int pivot[nterms];

      // QR decompose the model matrix
      for (i = 0; i < nterms; i++) pivot[i] = 0;
      pomp_qr(X,&nx,&nterms,pivot,tau);
      // solve R b = Q'y for b
      pomp_qrqy(yp,X,tau,&nx,&one,&nterms,&one,&one); // y <- Q'y 
      pomp_backsolve(X,&nx,&nterms,yp,b,&one); // b <- R^{-1} Q'y 
      
      // store b in kth column of beta
      for (i = 0; i < nterms; i++) beta[pivot[i]] = b[i]; 

    }

    Free(X);

  }

}
