// -*- C++ -*-

#include <Rdefines.h>
#include "pomp_internal.h"

// examines weights for filtering failure.
// computes conditional log likelihood and effective sample size.
// returns a named list.
// AT PRESENT, NO PARAMETER RESAMPLING IS PERFORMED.
SEXP wpfilter (SEXP X, SEXP Params, SEXP Weights, SEXP W, SEXP Trigger, SEXP Target, SEXP Np)
{

  SEXP dimX, dimP, Xnames, Pnames;
  int nvars, nreps, np;
  // REVISIT THIS IF PARAMETER RESAMPLING IS IMPLEMENTED
  //  int npars;
  //  int do_pr = 0;  // do parameter resampling?  at the moment, we do not
  int *dim, k;
  const char *dimnm[] = {"name",".id"};

  PROTECT(dimX = GET_DIM(X));
  dim = INTEGER(dimX);
  nvars = dim[0]; nreps = dim[1];
  PROTECT(Xnames = GET_ROWNAMES(GET_DIMNAMES(X)));

  PROTECT(Params = as_matrix(Params));
  PROTECT(dimP = GET_DIM(Params));
  dim = INTEGER(dimP);
  //  npars = dim[0];
  if (dim[1] > 1 && nreps != dim[1])
    err("ncol('params') does not match ncol('states')"); // # nocov
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(Params)));

  np = *INTEGER(AS_INTEGER(Np));
  // REVISIT THIS IF PARAMETER RESAMPLING IS IMPLEMENTED
  //  do_pr = (dim[1] > 1);	      // resample parameters as well as states

  PROTECT(Weights = duplicate(AS_NUMERIC(Weights)));
  PROTECT(W = duplicate(AS_NUMERIC(W)));

  int nprotect = 7;

  double trigger, target;
  trigger = *REAL(AS_NUMERIC(Trigger));
  target = *REAL(AS_NUMERIC(Target));

  // add W to Weights. check weights for NAN or +Inf. find max of weights.
  long double maxlogw = R_NegInf, oldw = 0;
  double *xW = REAL(Weights);
  double *xw = REAL(W);
  for (k = 0; k < nreps; k++, xw++, xW++) {

    oldw += exp(*xW);
    *xW += *xw;

    if (ISNAN(*xW) || *xW == R_PosInf) { // check the weights
      SEXP rv;
      PROTECT(rv = NEW_INTEGER(1)); nprotect++;
      *INTEGER(rv) = k+1; // return the index of the peccant particle
      UNPROTECT(nprotect);
      return rv;
    }

    if (maxlogw < *xW) maxlogw = *xW;

  }

  // set up return list
  SEXP retval, retvalnames, newdim;
  
  PROTECT(retval = NEW_LIST(5));
  PROTECT(retvalnames = NEW_CHARACTER(5));
  PROTECT(newdim = NEW_INTEGER(2));
  nprotect += 3;
    
  dim = INTEGER(newdim);
  dim[0] = nvars; dim[1] = nreps;
  SET_DIM(X,newdim);
  setrownames(X,Xnames,2);
  fixdimnames(X,dimnm,2);
  
  SET_STRING_ELT(retvalnames,0,mkChar("states"));
  SET_STRING_ELT(retvalnames,1,mkChar("params"));
  SET_STRING_ELT(retvalnames,2,mkChar("weights"));
  SET_STRING_ELT(retvalnames,3,mkChar("ess"));
  SET_STRING_ELT(retvalnames,4,mkChar("loglik"));
  SET_NAMES(retval,retvalnames);
  
  SEXP ess, loglik;
  PROTECT(ess = NEW_NUMERIC(1));    // effective sample size
  PROTECT(loglik = NEW_NUMERIC(1)); // conditional log likelihood
  nprotect += 2;

  if (maxlogw == R_NegInf) {	// all particles have zero likelihood
    *REAL(ess) = 0;
    *REAL(loglik) = R_NegInf;
  } else {
    // compute sum and sum of squares
    long double w = 0, ws = 0;
    xw = REAL(W);
    xW = REAL(Weights);
    for (k = 0; k < nreps; k++, xw++, xW++) {
      *xW -= maxlogw;
      *xw = exp(*xW);
      if (*xw != 0) {
        w += *xw;
        ws += (*xw)*(*xw);
      }
    }
    *(REAL(ess)) = w*w/ws;
    *(REAL(loglik)) = maxlogw + log(w) - log(oldw);
  }

  SET_ELEMENT(retval,3,ess);
  SET_ELEMENT(retval,4,loglik);
  
  if (maxlogw == R_NegInf || (*REAL(ess) >= nreps*trigger && np == nreps)) { // do not resample
    
    SET_ELEMENT(retval,0,X);
    SET_ELEMENT(retval,1,Params);
    SET_ELEMENT(retval,2,Weights);
    
  } else {                      // no resampling

    // create storage for new states
    SEXP newstates = R_NilValue;
    double *ss = 0, *st = 0;
    {
      int xdim[] = {nvars, np};
      PROTECT(newstates = makearray(2,xdim)); nprotect++;
      setrownames(newstates,Xnames,2);
      fixdimnames(newstates,dimnm,2);
      ss = REAL(X);
      st = REAL(newstates);
    }
    SET_ELEMENT(retval,0,newstates);
      
    // REVISIT THIS IF PARAMETER RESAMPLING IS IMPLEMENTED
    // create storage for new parameters
    //    SEXP newparams = R_NilValue;
    //    double *ps = 0, *pt = 0;
    // if (do_pr) {
    //   int xdim[] = {npars, np};
    //   PROTECT(newparams = makearray(2,xdim)); nprotect++;
    //   setrownames(newparams,Pnames,2);
    //   fixdimnames(newparams,dimnm,2);
    //   ps = REAL(Params);
    //   pt = REAL(newparams);
    //   SET_ELEMENT(retval,1,newparams);
    // } else {
    SET_ELEMENT(retval,1,Params);
    // }
      
    // create storage for new weights
    SEXP newweights = R_NilValue;
    double *ws = 0, *wt = 0;
    {
      PROTECT(newweights = NEW_NUMERIC(np)); nprotect++;
      ws = REAL(Weights);
      wt = REAL(newweights);
    }
    SET_ELEMENT(retval,2,newweights);
      
    // compute target resampling weights
    for (k = 0, xw = REAL(W), xW = REAL(Weights); k < nreps; k++, xw++, xW++) {
      *xw = *xW;
      *xW *= target;
      *xw = exp(*xw - *xW);
    }
    
    // do resampling
    {
      int sample[np];
      GetRNGstate();
      nosort_resamp(nreps,REAL(W),np,sample,0);
      PutRNGstate();
      
      for (k = 0; k < np; k++) { // copy the particles
	int sp = sample[k], j;
	double *xx;
	//	double *xp;
	for (j = 0, xx = ss+nvars*sp; j < nvars; j++, st++, xx++)
	  *st = *xx;
	// REVISIT THIS IF PARAMETER RESAMPLING IS IMPLEMENTED
	// if (do_pr) {
	//   for (j = 0, xp = ps+npars*sp; j < npars; j++, pt++, xp++)
	//     *pt = *xp;
	// }
	wt[k] = ws[sp];
      }
    }
  }

  UNPROTECT(nprotect);
  return(retval);
}
