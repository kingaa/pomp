// -*- C++ -*-

#include <Rdefines.h>
#include "pomp_internal.h"

// examines weights for filtering failure.
// computes conditional log likelihood and effective sample size.
// it is assumed that ncol(x) == ncol(params).
// returns a named list.
SEXP wpfilter (SEXP x, SEXP params, SEXP Np, SEXP weights)
{

  SEXP ess, loglik;
  SEXP retval, retvalnames;
  const char *dimnm[2] = {"variable","rep"};
  double *xw = 0;
  long double w = 0, ws, maxw;
  SEXP dimX, dimP, newdim, Xnames;
  int *dim;
  int nvars, nreps;
  int all_fail = 0;
  int k;

  PROTECT(dimX = GET_DIM(x));
  dim = INTEGER(dimX);
  nvars = dim[0]; nreps = dim[1];
  PROTECT(Xnames = GET_ROWNAMES(GET_DIMNAMES(x)));

  PROTECT(params = as_matrix(params));
  PROTECT(dimP = GET_DIM(params));
  dim = INTEGER(dimP);
  if (nreps % dim[1] != 0)
    err("ncol('states') should be a multiple of ncol('params')"); // # nocov

  PROTECT(weights = duplicate(weights)); // FIXME: unnecessary copy

  PROTECT(ess = NEW_NUMERIC(1));    // effective sample size
  PROTECT(loglik = NEW_NUMERIC(1)); // log likelihood

  int nprotect = 7;

  xw = REAL(weights);
  
  // check and normalize the weights
  for (k = 0, maxw = R_NegInf; k < nreps; k++) {

    if (ISNAN(xw[k]) || xw[k] == R_PosInf) { // check the weights
      PROTECT(retval = NEW_INTEGER(1)); nprotect++;
      *INTEGER(retval) = k+1; // return the index of the peccant particle
      UNPROTECT(nprotect);
      return retval;
    }    

    if (maxw < xw[k]) maxw = xw[k];

  }

  if (maxw == R_NegInf) all_fail = 1;

  if (all_fail) {
    
    *(REAL(loglik)) = R_NegInf;
    *(REAL(ess)) = 0;             // zero effective sample size

  } else {

    // compute sum and sum of squares
    for (k = 0, w = 0, ws = 0; k < nreps; k++) {
      xw[k] = exp(xw[k]-maxw);
      if (xw[k] != 0) {
	w += xw[k];
	ws += xw[k]*xw[k];
      }
    }

    *(REAL(loglik)) = maxw + log(w/((double) nreps)); // mean of weights is likelihood
    *(REAL(ess)) = w*w/ws;      // effective sample size
  }

  GetRNGstate();

  // don't resample: just drop 3rd dimension in x prior to return
  
  PROTECT(newdim = NEW_INTEGER(2)); nprotect++;
  dim = INTEGER(newdim);
  dim[0] = nvars; dim[1] = nreps;
  SET_DIM(x,newdim);
  setrownames(x,Xnames,2);
  fixdimnames(x,dimnm,2);

  PutRNGstate();

  PROTECT(retval = NEW_LIST(4));
  PROTECT(retvalnames = NEW_CHARACTER(4));
  nprotect += 2;
  SET_STRING_ELT(retvalnames,0,mkChar("loglik"));
  SET_STRING_ELT(retvalnames,1,mkChar("ess"));
  SET_STRING_ELT(retvalnames,2,mkChar("states"));
  SET_STRING_ELT(retvalnames,3,mkChar("params"));
  SET_NAMES(retval,retvalnames);

  SET_ELEMENT(retval,0,loglik);
  SET_ELEMENT(retval,1,ess);
  SET_ELEMENT(retval,2,x);
  SET_ELEMENT(retval,3,params);

  UNPROTECT(nprotect);
  return(retval);
}
