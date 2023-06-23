// -*- C++ -*-

#include "internal.h"
#include <Rdefines.h>

static void pred_mean_var (int, int, int, const double *, double *, double *);
static void filt_mean (int, int, int, long double, const double *, const double *, double *);

// computes conditional log likelihood and effective sample size.
// computes (if desired) prediction mean, prediction variance, filtering mean.
// it is assumed that ncol(x) == ncol(params).
// weights are used in filtering mean computation.
// if length(weights) == 1, an unweighted average is computed.
// tracks ancestry of particles if desired.
// returns all of the above in a named list.
SEXP pfilter (SEXP x, SEXP params, SEXP Np,
	      SEXP predmean, SEXP predvar,
	      SEXP filtmean, SEXP trackancestry, SEXP doparRS,
	      SEXP weights, SEXP wave)
{

  SEXP pm = R_NilValue, pv = R_NilValue, fm = R_NilValue;
  SEXP anc = R_NilValue, wmean = R_NilValue;
  SEXP ess, loglik;
  SEXP newstates = R_NilValue, newparams = R_NilValue;
  SEXP retval, retvalnames;
  const char *dimnm[2] = {"name",".id"};
  double *xw = 0;
  long double w = 0, ws, maxw, sum;
  int *xanc = 0;
  SEXP dimX, dimP, newdim, Xnames, Pnames;
  int *dim, np;
  int nvars, npars = 0, nreps;
  int do_pm, do_pv, do_fm, do_ta, do_pr, do_wave, all_fail = 0;
  int j, k;

  PROTECT(dimX = GET_DIM(x));
  dim = INTEGER(dimX);
  nvars = dim[0]; nreps = dim[1];
  PROTECT(Xnames = GET_ROWNAMES(GET_DIMNAMES(x)));

  PROTECT(params = as_matrix(params));
  PROTECT(dimP = GET_DIM(params));
  dim = INTEGER(dimP);
  npars = dim[0];
  if (nreps % dim[1] != 0)
    err("ncol('states') should be a multiple of ncol('params')"); // # nocov
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params)));

  PROTECT(weights = duplicate(weights));

  np = *(INTEGER(AS_INTEGER(Np))); // number of particles to resample

  do_pm = *(LOGICAL(AS_LOGICAL(predmean))); // calculate prediction means?
  do_pv = *(LOGICAL(AS_LOGICAL(predvar)));  // calculate prediction variances?
  do_fm = *(LOGICAL(AS_LOGICAL(filtmean))); // calculate filtering means?
  do_ta = *(LOGICAL(AS_LOGICAL(trackancestry))); // track ancestry?
  // Do we need to do parameter resampling?
  do_pr = *(LOGICAL(AS_LOGICAL(doparRS)));
  // Do we need to take a weighted average over the parameters?
  do_wave = *(LOGICAL(AS_LOGICAL(wave)));

  if (do_pr) {
    if (dim[1] != nreps)
      err("ncol('states') should be equal to ncol('params')"); // # nocov
  }

  PROTECT(ess = NEW_NUMERIC(1));    // effective sample size
  PROTECT(loglik = NEW_NUMERIC(1)); // log likelihood

  int nprotect = 8;

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

    *(REAL(loglik)) = maxw + log(w/((long double) nreps)); // mean of weights is likelihood
    *(REAL(ess)) = w*w/ws;      // effective sample size
  }

  if (do_pm || do_pv) {
    PROTECT(pm = NEW_NUMERIC(nvars)); nprotect++;
  }

  if (do_pv) {
    PROTECT(pv = NEW_NUMERIC(nvars)); nprotect++;
  }

  if (do_fm) {
    PROTECT(fm = NEW_NUMERIC(nvars)); nprotect++;
  }

  if (do_ta) {
    PROTECT(anc = NEW_INTEGER(np)); nprotect++;
    xanc = INTEGER(anc);
  }

  if (do_wave) {
    PROTECT(wmean = NEW_NUMERIC(npars)); nprotect++;
    SET_NAMES(wmean,Pnames);
  }

  if (do_pm || do_pv) {
    double *tmp = (do_pv) ? REAL(pv) : 0;
    pred_mean_var(nvars,nreps,do_pv,REAL(x),REAL(pm),tmp);
  }
  
  //  compute filter mean
  if (do_fm) {
    filt_mean(nvars,nreps,all_fail,w,xw,REAL(x),REAL(fm));
  }

  // compute weighted average of parameters
  if (do_wave) {
    if (all_fail)
      warn("%s %s","filtering failure at last filter iteration:",
	   "using unweighted mean for point estimate.");
    double *xwm = REAL(wmean);
    for (j = 0; j < npars; j++, xwm++) {
      double *xp = REAL(params)+j;
      if (all_fail) {           // unweighted average
        for (k = 0, sum = 0; k < nreps; k++, xp += npars) sum += *xp;
        *xwm = sum/((long double) nreps);
      } else {
        for (k = 0, sum = 0; k < nreps; k++, xp += npars) {
	  if (xw[k]!=0) sum += xw[k]*(*xp);
	}
        *xwm = sum/w;
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
      int sp = sample[k];
      double *xx, *xp;
      for (j = 0, xx = ss+nvars*sp; j < nvars; j++, st++, xx++)
        *st = *xx;
      if (do_pr) {
        for (j = 0, xp = ps+npars*sp; j < npars; j++, pt++, xp++)
          *pt = *xp;
      }
      if (do_ta) xanc[k] = sp+1;
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

  PROTECT(retval = NEW_LIST(9));
  PROTECT(retvalnames = NEW_CHARACTER(9));
  nprotect += 2;
  SET_STRING_ELT(retvalnames,0,mkChar("loglik"));
  SET_STRING_ELT(retvalnames,1,mkChar("ess"));
  SET_STRING_ELT(retvalnames,2,mkChar("states"));
  SET_STRING_ELT(retvalnames,3,mkChar("params"));
  SET_STRING_ELT(retvalnames,4,mkChar("pm"));
  SET_STRING_ELT(retvalnames,5,mkChar("pv"));
  SET_STRING_ELT(retvalnames,6,mkChar("fm"));
  SET_STRING_ELT(retvalnames,7,mkChar("ancestry"));
  SET_STRING_ELT(retvalnames,8,mkChar("wmean"));
  SET_NAMES(retval,retvalnames);

  SET_ELEMENT(retval,0,loglik);
  SET_ELEMENT(retval,1,ess);

  if (all_fail) {
    SET_ELEMENT(retval,2,x);
  } else {
    SET_ELEMENT(retval,2,newstates);
  }

  if (all_fail || !do_pr) {
    SET_ELEMENT(retval,3,params);
  } else {
    SET_ELEMENT(retval,3,newparams);
  }

  if (do_pm) {
    SET_ELEMENT(retval,4,pm);
  }
  if (do_pv) {
    SET_ELEMENT(retval,5,pv);
  }
  if (do_fm) {
    SET_ELEMENT(retval,6,fm);
  }
  if (do_ta) {
    SET_ELEMENT(retval,7,anc);
  }
  if (do_wave) {
    SET_ELEMENT(retval,8,wmean);
  }

  UNPROTECT(nprotect);
  return(retval);
}

static void pred_mean_var (int nvars, int nreps, int do_pv,
			   const double *x,
			   double *pm, double *pv)
{
  long double sum, sumsq, vsq;
  const double *xx;
  int j, k;

  for (j = 0; j < nvars; j++, pm++, pv++) {

    xx = x+j;
    
    // compute prediction mean
    for (k = 0, sum = 0; k < nreps; k++, xx += nvars) sum += *xx;
    *pm = sum/((long double) nreps);
    
    // compute prediction variance
    if (do_pv) {
      xx = x+j;
      for (k = 0, sumsq = 0; k < nreps; k++, xx += nvars) {
	vsq = *xx - sum;
	sumsq += vsq*vsq;
      }
      *pv = sumsq/((long double) (nreps - 1));
    }
  }
}

static void filt_mean (int nvars, int nreps, int all_fail, long double wsum,
		       const double *w, const double *x, double *fm)
{
  long double sum;
  const double *xx;
  int j, k;
  for (j = 0; j < nvars; j++, fm++) {
    xx = x+j;
    if (all_fail) {           // unweighted average
      for (k = 0, sum = 0; k < nreps; k++, xx += nvars) sum += *xx;
      *fm = sum/((long double) nreps);
    } else {                  // weighted average
      for (k = 0, sum = 0; k < nreps; k++, xx += nvars) sum += w[k]*(*xx);
      *fm = sum/wsum;
    }
  }
}
