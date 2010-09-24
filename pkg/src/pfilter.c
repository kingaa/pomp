// -*- C++ -*-

#include "pomp_internal.h"
#include <Rdefines.h>

// examines weights for filtering failure
// computes log likelihood and effective sample size
// computes (if desired) prediction mean, prediction variance, filtering mean.
// it is assumed that ncol(x) == ncol(params).
// weights are used in filtering mean computation.
// if length(weights) == 1, an unweighted average is computed.
// returns all of the above in a named list
SEXP pfilter_computations (SEXP x, SEXP params, 
			   SEXP rw, SEXP rw_names, 
			   SEXP predmean, SEXP predvar,
			   SEXP filtmean, SEXP weights, SEXP tol)
{
  int nprotect = 0;
  SEXP pm, pv, fm, ess, fail, loglik;
  SEXP retval, retvalnames;
  double *xpm, *xpv, *xfm, *xw;
  SEXP dimX, Pnames, pindex;
  double *xx, *xp;
  int *dim, *pidx, lv;
  int nvars, npars, nrw, nreps, offset, nlost;
  int do_rw, do_pm, do_pv, do_fm, all_fail = 0;
  double sum, sumsq, vsq, ws, w, toler;
  int j, k;

  PROTECT(dimX = GET_DIM(x)); nprotect++;
  dim = INTEGER(dimX);
  nvars = dim[0]; nreps = dim[1];
  xx = REAL(x);
  do_rw = *(LOGICAL(AS_LOGICAL(rw)));
  do_pm = *(LOGICAL(AS_LOGICAL(predmean)));
  do_pv = *(LOGICAL(AS_LOGICAL(predvar)));
  do_fm = *(LOGICAL(AS_LOGICAL(filtmean)));

  PROTECT(ess = NEW_NUMERIC(1)); nprotect++;
  PROTECT(loglik = NEW_NUMERIC(1)); nprotect++;
  PROTECT(fail = NEW_LOGICAL(1)); nprotect++;

  xw = REAL(weights);
  toler = *(REAL(tol));

  // check the weights and compute sum and sum of squares
  for (k = 0, w = 0, ws = 0, nlost = 0; k < nreps; k++) {
    if (xw[k] > toler) {	
      w += xw[k];
      ws += xw[k]*xw[k];
    } else {			// this particle is lost
      xw[k] = 0;
      nlost++;
    }
  }
  if (nlost >= nreps) all_fail = 1; // all particles are lost
  if (all_fail) {
    *(REAL(loglik)) = log(toler); // minimum log-likelihood
    *(REAL(ess)) = 0;		    // zero effective sample size
  } else {
    *(REAL(loglik)) = log(w/((double) nreps)); // mean of weights is likelihood
    *(REAL(ess)) = w*w/ws;			 // effective sample size
  }
  *(LOGICAL(fail)) = all_fail;

  if (do_rw) {
    PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;
    PROTECT(pindex = matchnames(Pnames,rw_names)); nprotect++;
    xp = REAL(params);
    pidx = INTEGER(pindex);
    npars = LENGTH(Pnames);
    nrw = LENGTH(rw_names);
    lv = nvars+nrw;
  } else {
    pidx = NULL;
    lv = nvars;
  }

  if (do_pm || do_pv) {
    PROTECT(pm = NEW_NUMERIC(lv)); nprotect++;
    xpm = REAL(pm);
  }

  if (do_pv) {
    PROTECT(pv = NEW_NUMERIC(lv)); nprotect++;
    xpv = REAL(pv);
  }

  if (do_fm) {
    if (do_rw) {
      PROTECT(fm = NEW_NUMERIC(nvars+npars)); nprotect++;
    } else {
      PROTECT(fm = NEW_NUMERIC(nvars)); nprotect++;
    }
    xfm = REAL(fm);
  }

  for (j = 0; j < nvars; j++) {	// state variables

    // compute prediction mean
    if (do_pm || do_pv) {
      for (k = 0, sum = 0; k < nreps; k++) sum += xx[j+k*nvars];
      sum /= ((double) nreps);
      xpm[j] = sum;
    }

    // compute prediction variance
    if (do_pv) {	
      for (k = 0, sumsq = 0; k < nreps; k++) {
	vsq = xx[j+k*nvars]-sum;
	sumsq += vsq*vsq;
      }
      xpv[j] = sumsq / ((double) (nreps - 1));
    }

    //  compute filter mean
    if (do_fm) {
      if (all_fail) {
	for (k = 0, ws = 0; k < nreps; k++) ws += xx[j+k*nvars]; 
	xfm[j] = ws/((double) nreps);
      } else { 
	for (k = 0, ws = 0; k < nreps; k++) ws += xx[j+k*nvars]*xw[k]; 
	xfm[j] = ws/w;
      }
    }

  }

  // compute means and variances for parameters
  if (do_rw) {
    for (j = 0; j < nrw; j++) {
      offset = pidx[j];

      if (do_pm || do_pv) {
	for (k = 0, sum = 0; k < nreps; k++) sum += xp[offset+k*npars];
	sum /= ((double) nreps);
	xpm[nvars+j] = sum;
      }

      if (do_pv) {
	for (k = 0, sumsq = 0; k < nreps; k++) {
	  vsq = xp[offset+k*npars]-sum;
	  sumsq += vsq*vsq;
	}
	xpv[nvars+j] = sumsq / ((double) (nreps - 1));
      }

    }
    if (do_fm) {
      for (j = 0; j < npars; j++) {
	if (all_fail) {
	  for (k = 0, ws = 0; k < nreps; k++) ws += xp[j+k*npars];
	  xfm[nvars+j] = ws/((double) nreps);
	} else {
	  for (k = 0, ws = 0; k < nreps; k++) ws += xp[j+k*npars]*xw[k];
	  xfm[nvars+j] = ws/w;
	}
      }
    }
  }

  PROTECT(retval = NEW_LIST(6)); nprotect++;
  PROTECT(retvalnames = NEW_CHARACTER(6)); nprotect++;
  SET_STRING_ELT(retvalnames,0,mkChar("fail"));
  SET_STRING_ELT(retvalnames,1,mkChar("loglik"));
  SET_STRING_ELT(retvalnames,2,mkChar("ess"));
  SET_STRING_ELT(retvalnames,3,mkChar("pm"));
  SET_STRING_ELT(retvalnames,4,mkChar("pv"));
  SET_STRING_ELT(retvalnames,5,mkChar("fm"));
  SET_NAMES(retval,retvalnames);

  SET_ELEMENT(retval,0,fail);
  SET_ELEMENT(retval,1,loglik);
  SET_ELEMENT(retval,2,ess);

  if (do_pm) {
    SET_ELEMENT(retval,3,pm);
  }
  if (do_pv) {
    SET_ELEMENT(retval,4,pv);
  }
  if (do_fm) {
    SET_ELEMENT(retval,5,fm);
  }

  UNPROTECT(nprotect);
  return(retval);
}

