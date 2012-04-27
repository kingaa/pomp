// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "pomp_internal.h"

SEXP do_skeleton (SEXP object, SEXP x, SEXP t, SEXP params, SEXP fun)
{
  int nprotect = 0;
  int nvars, npars, nrepp, nrepx, nreps, ntimes, ncovars;
  int mode;
  int use_names;
  int *dim, ndim[3], *op;
  SEXP fn, fcall, rho, ans, nm;
  SEXP tvec, xvec, pvec, cvec;
  SEXP Snames, Cnames, Pnames;
  SEXP statenames, paramnames, covarnames;
  SEXP F;
  double *xs, *ts, *ps, *fs, *ft;
  double *tp, *xp, *pp, *cp;
  int *sidx, *pidx, *cidx;
  int first = 1;
  pomp_skeleton *ff = NULL;
  struct lookup_table covariate_table;
  int i, j, k;

  PROTECT(t = AS_NUMERIC(t)); nprotect++;
  ntimes = LENGTH(t);
  ts = REAL(t);

  PROTECT(x = as_state_array(x)); nprotect++;
  dim = INTEGER(GET_DIM(x));
  nvars = dim[0]; nrepx = dim[1];
  if (ntimes != dim[2])
    error("skeleton error: length of 't' and 3rd dimension of 'x' do not agree");
  xs = REAL(x);

  PROTECT(params = as_matrix(params)); nprotect++;
  dim = INTEGER(GET_DIM(params));
  npars = dim[0]; nrepp = dim[1];
  ps = REAL(params);

  // 2nd dimension of 'x' and 'params' need not agree
  nreps = (nrepp > nrepx) ? nrepp : nrepx;
  if ((nreps % nrepp != 0) || (nreps % nrepx != 0))
    error("skeleton error: 2nd dimensions of 'x' and 'params' are incompatible");

  PROTECT(Snames = GET_ROWNAMES(GET_DIMNAMES(x))); nprotect++;
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;
  PROTECT(Cnames = GET_COLNAMES(GET_DIMNAMES(GET_SLOT(object,install("covar"))))); nprotect++;
    
  // set up the covariate table
  covariate_table = make_covariate_table(object,&ncovars);

  // vector for interpolated covariates
  PROTECT(cvec = NEW_NUMERIC(ncovars)); nprotect++;
  SET_NAMES(cvec,Cnames);
  cp = REAL(cvec);

  // set up the array to be returned
  ndim[0] = nvars; ndim[1] = nreps; ndim[2] = ntimes;
  PROTECT(F = makearray(3,ndim)); nprotect++; 
  setrownames(F,Snames,3);

  // extract the user-defined function
  PROTECT(fn = unpack_pomp_fun(fun,&mode)); nprotect++;

  // extract 'userdata' as pairlist
  PROTECT(fcall = VectorToPairList(GET_SLOT(object,install("userdata")))); nprotect++;

  // first do setup
  switch (mode) {
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
    
    // construct state, parameter, covariate, observable indices
    sidx = INTEGER(PROTECT(name_index(Snames,object,"statenames"))); nprotect++;
    pidx = INTEGER(PROTECT(name_index(Pnames,object,"paramnames"))); nprotect++;
    cidx = INTEGER(PROTECT(name_index(Cnames,object,"covarnames"))); nprotect++;

    ff = (pomp_skeleton *) R_ExternalPtrAddr(fn);
    
    break;

  default:
    error("unrecognized 'mode' slot in 'skeleton'");
    break;
  }


  // now do computations
  switch (mode) {
  case 0: 			// R skeleton

    for (k = 0, ft = REAL(F); k < ntimes; k++, ts++) { // loop over times

      R_CheckUserInterrupt();	// check for user interrupt
      
      *tp = *ts;		// copy the time

      // interpolate the covar functions for the covariates
      table_lookup(&covariate_table,*tp,cp,0);
      
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
    
    set_pomp_userdata(fcall);

    for (k = 0, ft = REAL(F); k < ntimes; k++, ts++) { // loop over times

      R_CheckUserInterrupt();	// check for user interrupt
      
      // interpolate the covar functions for the covariates
      table_lookup(&covariate_table,*ts,cp,0);
      
      for (j = 0; j < nreps; j++, ft += nvars) { // loop over replicates
	
	xp = &xs[nvars*((j%nrepx)+nrepx*k)];
	pp = &ps[npars*(j%nrepp)];
	
	(*ff)(ft,xp,pp,sidx,pidx,cidx,ncovars,cp,*ts);
	
      }
    }

    unset_pomp_userdata();

    break;

  default:
    error("unrecognized 'mode' slot in 'skeleton'");
    break;
  }

  UNPROTECT(nprotect);
  return F;
}
