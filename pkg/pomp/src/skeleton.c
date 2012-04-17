// dear emacs, please treat this as -*- C++ -*-

#include "pomp_internal.h"
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP do_skeleton (SEXP object, SEXP x, SEXP t, SEXP params, SEXP fun)
{
  int nprotect = 0;
  int nvars, npars, nrepp, nrepx, nreps, ntimes, covlen, covdim;
  int use_native;
  int fdim[3];
  SEXP tcovar, covar, fn, Snames, F;
  double *xs, *ts, *ps, *ft;
  int *dim;

  SEXP rho, fcall, ans, nm;
  SEXP tvec, xvec, pvec, cvec;
  int first = 1;
  int use_names = 0;
  int *op;

  int *sidx, *pidx, *cidx;
  SEXP statenames, paramnames, covarnames;
  pomp_skeleton *vf;

  double *tp, *xp, *pp, *cp, *fs;
  int i, j, k;

  PROTECT(t = AS_NUMERIC(t)); nprotect++;
  ntimes = LENGTH(t);
  ts = REAL(t);

  PROTECT(x = as_state_array(x)); nprotect++;
  dim = INTEGER(GET_DIM(x));
  nvars = dim[0]; nrepx = dim[1];
  if (ntimes != dim[2])
    error("skeleton error: length of 't' and 3rd dimension of 'x' do not agree");
  PROTECT(Snames = GET_ROWNAMES(GET_DIMNAMES(x))); nprotect++;
  xs = REAL(x);

  PROTECT(params = as_matrix(params)); nprotect++;
  dim = INTEGER(GET_DIM(params));
  npars = dim[0]; nrepp = dim[1];
  ps = REAL(params);

  // 2nd dimension of 'x' and 'params' need not agree
  nreps = (nrepp > nrepx) ? nrepp : nrepx;
  if ((nreps % nrepp != 0) || (nreps % nrepx != 0))
    error("skeleton error: 2nd dimensions of 'x' and 'params' are incompatible");

  // extract the covariates
  PROTECT(tcovar =  GET_SLOT(object,install("tcovar"))); nprotect++;
  PROTECT(covar =  GET_SLOT(object,install("covar"))); nprotect++;
  dim = INTEGER(GET_DIM(covar)); 
  covlen = dim[0]; covdim = dim[1];

  // set up the covariate table
  struct lookup_table covariate_table = {covlen, covdim, 0, REAL(tcovar), REAL(covar)};

  // set up the array to be returned
  fdim[0] = nvars; fdim[1] = nreps; fdim[2] = ntimes;
  PROTECT(F = makearray(3,fdim)); nprotect++; 
  setrownames(F,Snames,3);

  // extract the user-defined function
  PROTECT(fn = VECTOR_ELT(fun,0)); nprotect++;
  use_native = INTEGER(VECTOR_ELT(fun,1))[0];
  
  PROTECT(cvec = NEW_NUMERIC(covdim)); nprotect++;
  SET_NAMES(cvec,GET_COLNAMES(GET_DIMNAMES(covar)));
  cp = REAL(cvec);

  // first do setup
  switch (use_native) {
  case 0: 			// R skeleton

    PROTECT(tvec = NEW_NUMERIC(1)); nprotect++;
    tp = REAL(tvec);

    PROTECT(xvec = NEW_NUMERIC(nvars)); nprotect++;
    SET_NAMES(xvec,Snames);
    xp = REAL(xvec);

    PROTECT(pvec = NEW_NUMERIC(npars)); nprotect++;
    SET_NAMES(pvec,GET_ROWNAMES(GET_DIMNAMES(params)));
    pp = REAL(pvec);

    // set up the function call
    PROTECT(fcall = VectorToPairList(GET_SLOT(object,install("userdata")))); nprotect++;
    PROTECT(fcall = LCONS(cvec,fcall)); nprotect++;
    SET_TAG(fcall,install("covars"));
    PROTECT(fcall = LCONS(pvec,fcall)); nprotect++;
    SET_TAG(fcall,install("params"));
    PROTECT(fcall = LCONS(tvec,fcall)); nprotect++;
    SET_TAG(fcall,install("t"));
    PROTECT(fcall = LCONS(xvec,fcall)); nprotect++;
    SET_TAG(fcall,install("x"));
    PROTECT(fcall = LCONS(fn,fcall)); nprotect++;
    PROTECT(rho = (CLOENV(fn))); nprotect++;

    break;

  case 1:			// native skeleton
    
    vf = (pomp_skeleton *) R_ExternalPtrAddr(fn);
    
    PROTECT(statenames = GET_SLOT(object,install("statenames"))); nprotect++;
    if (LENGTH(statenames) > 0) {
      sidx = INTEGER(PROTECT(MATCHROWNAMES(x,statenames))); nprotect++;
    } else {
      sidx = 0;
    }
    
    PROTECT(paramnames = GET_SLOT(object,install("paramnames"))); nprotect++;
    if (LENGTH(paramnames) > 0) {
      pidx = INTEGER(PROTECT(MATCHROWNAMES(params,paramnames))); nprotect++;
    } else {
      pidx = 0;
    }
    
    PROTECT(covarnames = GET_SLOT(object,install("covarnames"))); nprotect++;
    if (LENGTH(covarnames) > 0) {
      cidx = INTEGER(PROTECT(MATCHCOLNAMES(covar,covarnames))); nprotect++;
    } else {
      cidx = 0;
    }

    break;

  default:
    error("unrecognized 'use' slot in 'skeleton'");
    break;
  }


  // now do computations
  switch (use_native) {
  case 0: 			// R skeleton

    for (k = 0, ft = REAL(F); k < ntimes; k++, ts++) { // loop over times

      R_CheckUserInterrupt();	// check for user interrupt
      
      *tp = *ts;		// copy the time

      // interpolate the covar functions for the covariates
      if (covdim > 0) table_lookup(&covariate_table,*ts,cp,0);
      
      for (j = 0; j < nreps; j++, ft += nvars) { // loop over replicates
	
	for (i = 0; i < nvars; i++) xp[i] = xs[i+nvars*((j%nrepx)+nrepx*k)];
	for (i = 0; i < npars; i++) pp[i] = ps[i+npars*(j%nrepp)];
	
	if (first) {
	  
	  PROTECT(ans = eval(fcall,rho)); nprotect++;
	  if (LENGTH(ans)!=nvars)
	    error("user 'skeleton' returns a vector of %d state variables but %d are expected",LENGTH(ans),nvars);

	  // get name information to fix possible alignment problems
	  PROTECT(nm = GET_NAMES(ans)); nprotect++;
	  use_names = !isNull(nm);
	  if (use_names) {
	    op = INTEGER(PROTECT(matchnames(Snames,nm))); nprotect++;
	  } else {
	    op = 0;
	  }
	  
	  fs = REAL(AS_NUMERIC(ans));
	  
	  first = 0;
	  
	} else {
	  
	  fs = REAL(AS_NUMERIC(eval(fcall,rho)));

	}
	  
	if (use_names) 
	  for (i = 0; i < nvars; i++) ft[op[i]] = fs[i];
	else
	  for (i = 0; i < nvars; i++) ft[i] = fs[i];
	
      }
    }

    break;

  case 1:			// native skeleton
    
    for (k = 0, ft = REAL(F); k < ntimes; k++, ts++) { // loop over times

      R_CheckUserInterrupt();	// check for user interrupt
      
      // interpolate the covar functions for the covariates
      if (covdim > 0) table_lookup(&covariate_table,*ts,cp,0);
      
      for (j = 0; j < nreps; j++, ft += nvars) { // loop over replicates
	
	xp = &xs[nvars*((j%nrepx)+nrepx*k)];
	pp = &ps[npars*(j%nrepp)];
	
	(*vf)(ft,xp,pp,sidx,pidx,cidx,covdim,cp,*ts);
	
      }
    }

    break;

  default:
    error("unrecognized 'use' slot in 'skeleton'");
    break;
  }

  UNPROTECT(nprotect);
  return F;
}
