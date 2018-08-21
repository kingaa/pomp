// -*- C++ -*-

#include "pomp_internal.h"
#include <Rdefines.h>

// examines weights for filtering failure.
// computes conditional log likelihood and effective sample size.
// computes (if desired) prediction mean, prediction variance, filtering mean.
// it is assumed that ncol(x) == ncol(params).
// weights are used in filtering mean computation.
// if length(weights) == 1, an unweighted average is computed.
// tracks ancestry of particles if desired.
// returns all of the above in a named list.
SEXP pfilter_computations (SEXP x, SEXP params, SEXP Np,
  SEXP predmean, SEXP predvar,
  SEXP filtmean, SEXP trackancestry, SEXP doparRS,
  SEXP weights, SEXP tol)
{
  int nprotect = 0;
  SEXP pm = R_NilValue, pv = R_NilValue, fm = R_NilValue, anc = R_NilValue;
  SEXP ess, fail, loglik;
  SEXP newstates = R_NilValue, newparams = R_NilValue;
  SEXP retval, retvalnames;
  const char *dimnm[2] = {"variable","rep"};
  double *xpm = 0, *xpv = 0, *xfm = 0, *xw = 0, *xx = 0, *xp = 0;
  int *xanc = 0;
  SEXP dimX, dimP, newdim, Xnames, Pnames;
  int *dim, np;
  int nvars, npars = 0, nreps, nlost;
  int do_pm, do_pv, do_fm, do_ta, do_pr, all_fail = 0;
  double sum = 0, sumsq = 0, vsq, ws, w, toler;
  int j, k;

  PROTECT(dimX = GET_DIM(x)); nprotect++;
  dim = INTEGER(dimX);
  nvars = dim[0]; nreps = dim[1];
  xx = REAL(x);
  PROTECT(Xnames = GET_ROWNAMES(GET_DIMNAMES(x))); nprotect++;

  PROTECT(params = as_matrix(params)); nprotect++;
  PROTECT(dimP = GET_DIM(params)); nprotect++;
  dim = INTEGER(dimP);
  npars = dim[0];
  if (nreps % dim[1] != 0)
    errorcall(R_NilValue,"ncol('states') should be a multiple of ncol('params')"); // # nocov
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;

  np = *(INTEGER(AS_INTEGER(Np))); // number of particles to resample

  do_pm = *(LOGICAL(AS_LOGICAL(predmean))); // calculate prediction means?
  do_pv = *(LOGICAL(AS_LOGICAL(predvar)));  // calculate prediction variances?
  do_fm = *(LOGICAL(AS_LOGICAL(filtmean))); // calculate filtering means?
  do_ta = *(LOGICAL(AS_LOGICAL(trackancestry))); // track ancestry?
  // Do we need to do parameter resampling?
  do_pr = *(LOGICAL(AS_LOGICAL(doparRS)));

  if (do_pr) {
    if (dim[1] != nreps)
      errorcall(R_NilValue,"ncol('states') should be equal to ncol('params')"); // # nocov
  }

  PROTECT(ess = NEW_NUMERIC(1)); nprotect++; // effective sample size
  PROTECT(loglik = NEW_NUMERIC(1)); nprotect++; // log likelihood
  PROTECT(fail = NEW_LOGICAL(1)); nprotect++;	// particle failure?

  xw = REAL(weights);
  toler = *(REAL(tol));		// failure tolerance

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
    *(REAL(ess)) = 0;		  // zero effective sample size
  } else {
    *(REAL(loglik)) = log(w/((double) nreps)); // mean of weights is likelihood
    *(REAL(ess)) = w*w/ws;	// effective sample size
  }
  *(LOGICAL(fail)) = all_fail;

  if (do_pm || do_pv) {
    PROTECT(pm = NEW_NUMERIC(nvars)); nprotect++;
    xpm = REAL(pm);
  }

  if (do_pv) {
    PROTECT(pv = NEW_NUMERIC(nvars)); nprotect++;
    xpv = REAL(pv);
  }

  if (do_fm) {
    PROTECT(fm = NEW_NUMERIC(nvars)); nprotect++;
    xfm = REAL(fm);
  }

  if (do_ta) {
    PROTECT(anc = NEW_INTEGER(np)); nprotect++;
    xanc = INTEGER(anc);
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
      if (all_fail) {		// unweighted average
        for (k = 0, ws = 0; k < nreps; k++) ws += xx[j+k*nvars];
        xfm[j] = ws/((double) nreps);
      } else { 			// weighted average
        for (k = 0, ws = 0; k < nreps; k++) ws += xx[j+k*nvars]*xw[k];
        xfm[j] = ws/w;
      }
    }

  }

  GetRNGstate();

  if (!all_fail) { // resample the particles unless we have filtering failure
    int xdim[2];
    int sample[np];
    double *ss = 0, *st = 0, *ps = 0, *pt = 0;

    // create storage for new states
    xdim[0] = nvars; xdim[1] = np;
    PROTECT(newstates = makearray(2,xdim)); nprotect++;
    setrownames(newstates,Xnames,2);
    fixdimnames(newstates,dimnm,2);
    ss = REAL(x);
    st = REAL(newstates);

    // create storage for new parameters
    if (do_pr) {
      xdim[0] = npars; xdim[1] = np;
      PROTECT(newparams = makearray(2,xdim)); nprotect++;
      setrownames(newparams,Pnames,2);
      fixdimnames(newparams,dimnm,2);
      ps = REAL(params);
      pt = REAL(newparams);
    }

    // resample
    nosort_resamp(nreps,REAL(weights),np,sample,0);
    for (k = 0; k < np; k++) { // copy the particles
      for (j = 0, xx = ss+nvars*sample[k]; j < nvars; j++, st++, xx++)
        *st = *xx;
      if (do_pr) {
        for (j = 0, xp = ps+npars*sample[k]; j < npars; j++, pt++, xp++)
          *pt = *xp;
      }
      if (do_ta) xanc[k] = sample[k]+1;
    }

  } else { // don't resample: just drop 3rd dimension in x prior to return

    PROTECT(newdim = NEW_INTEGER(2)); nprotect++;
    dim = INTEGER(newdim);
    dim[0] = nvars; dim[1] = nreps;
    SET_DIM(x,newdim);
    setrownames(x,Xnames,2);
    fixdimnames(x,dimnm,2);

    if (do_ta)
      for (k = 0; k < np; k++) xanc[k] = k+1;
  }

  PutRNGstate();

  PROTECT(retval = NEW_LIST(9)); nprotect++;
  PROTECT(retvalnames = NEW_CHARACTER(9)); nprotect++;
  SET_STRING_ELT(retvalnames,0,mkChar("fail"));
  SET_STRING_ELT(retvalnames,1,mkChar("loglik"));
  SET_STRING_ELT(retvalnames,2,mkChar("ess"));
  SET_STRING_ELT(retvalnames,3,mkChar("states"));
  SET_STRING_ELT(retvalnames,4,mkChar("params"));
  SET_STRING_ELT(retvalnames,5,mkChar("pm"));
  SET_STRING_ELT(retvalnames,6,mkChar("pv"));
  SET_STRING_ELT(retvalnames,7,mkChar("fm"));
  SET_STRING_ELT(retvalnames,8,mkChar("ancestry"));
  SET_NAMES(retval,retvalnames);

  SET_ELEMENT(retval,0,fail);
  SET_ELEMENT(retval,1,loglik);
  SET_ELEMENT(retval,2,ess);

  if (all_fail) {
    SET_ELEMENT(retval,3,x);
  } else {
    SET_ELEMENT(retval,3,newstates);
  }

  if (all_fail || !do_pr) {
    SET_ELEMENT(retval,4,params);
  } else {
    SET_ELEMENT(retval,4,newparams);
  }

  if (do_pm) {
    SET_ELEMENT(retval,5,pm);
  }
  if (do_pv) {
    SET_ELEMENT(retval,6,pv);
  }
  if (do_fm) {
    SET_ELEMENT(retval,7,fm);
  }
  if (do_ta) {
    SET_ELEMENT(retval,8,anc);
  }

  UNPROTECT(nprotect);
  return(retval);
}
