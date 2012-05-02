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
  int mode = -1;
  int *dim;
  SEXP fn, fcall, rho, ans, nm;
  SEXP tvec, xvec, pvec, cvec;
  SEXP Snames, Cnames, Pnames;
  SEXP F;
  int *sidx = 0, *pidx = 0, *cidx = 0;
  pomp_skeleton *ff = NULL;
  struct lookup_table covariate_table;

  PROTECT(t = AS_NUMERIC(t)); nprotect++;
  ntimes = LENGTH(t);

  PROTECT(x = as_state_array(x)); nprotect++;
  dim = INTEGER(GET_DIM(x));
  nvars = dim[0]; nrepx = dim[1];
  if (ntimes != dim[2])
    error("skeleton error: length of 't' and 3rd dimension of 'x' do not agree");

  PROTECT(params = as_matrix(params)); nprotect++;
  dim = INTEGER(GET_DIM(params));
  npars = dim[0]; nrepp = dim[1];

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

  // extract the user-defined function
  PROTECT(fn = unpack_pomp_fun(fun,&mode)); nprotect++;

  // extract 'userdata' as pairlist
  PROTECT(fcall = VectorToPairList(GET_SLOT(object,install("userdata")))); nprotect++;

  // first do setup
  switch (mode) {
  case 0: 			// R skeleton

    PROTECT(tvec = NEW_NUMERIC(1)); nprotect++;
    PROTECT(xvec = NEW_NUMERIC(nvars)); nprotect++;
    PROTECT(pvec = NEW_NUMERIC(npars)); nprotect++;
    SET_NAMES(xvec,Snames);
    SET_NAMES(pvec,Pnames);

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


  // set up the array to hold results
  {
    int dim[3] = {nvars, nreps, ntimes};
    PROTECT(F = makearray(3,dim)); nprotect++; 
    setrownames(F,Snames,3);
  }

  // now do computations
  switch (mode) {
  case 0: 			// R skeleton

    {
      int first = 1;
      int use_names;
      double *time = REAL(t);
      double *xs = REAL(x);
      double *ps = REAL(params);
      double *cp = REAL(cvec);
      double *tp = REAL(tvec);
      double *xp = REAL(xvec);
      double *pp = REAL(pvec);
      double *ft = REAL(F);
      double *fs;
      int *posn;
      int i, j, k;

      for (k = 0; k < ntimes; k++, time++) { // loop over times

	R_CheckUserInterrupt();	// check for user interrupt
      
	*tp = *time;		// copy the time

	// interpolate the covar functions for the covariates
	table_lookup(&covariate_table,*time,cp,0);
      
	for (j = 0; j < nreps; j++, ft += nvars) { // loop over replicates
	
	  memcpy(xp,&xs[nvars*((j%nrepx)+nrepx*k)],nvars*sizeof(double));
	  memcpy(pp,&ps[npars*(j%nrepp)],npars*sizeof(double));
	
	  if (first) {
	  
	    PROTECT(ans = eval(fcall,rho)); nprotect++;
	    if (LENGTH(ans)!=nvars)
	      error("user 'skeleton' returns a vector of %d state variables but %d are expected",LENGTH(ans),nvars);

	    // get name information to fix possible alignment problems
	    PROTECT(nm = GET_NAMES(ans)); nprotect++;
	    use_names = !isNull(nm);
	    if (use_names) {
	      posn = INTEGER(PROTECT(matchnames(Snames,nm))); nprotect++;
	    } else {
	      posn = 0;
	    }
	  
	    fs = REAL(AS_NUMERIC(ans));
	  
	    first = 0;
	  
	  } else {
	  
	    fs = REAL(AS_NUMERIC(eval(fcall,rho)));

	  }
	  
	  if (use_names) 
	    for (i = 0; i < nvars; i++) ft[posn[i]] = fs[i];
	  else
	    for (i = 0; i < nvars; i++) ft[i] = fs[i];
	
	}
      }
    }

    break;

  case 1:			// native skeleton
    
    set_pomp_userdata(fcall);

    {
      double *time = REAL(t);
      double *xs = REAL(x);
      double *ps = REAL(params);
      double *cp = REAL(cvec);
      double *ft = REAL(F);
      double *xp, *pp;
      int j, k;

      for (k = 0; k < ntimes; k++, time++) { // loop over times

	R_CheckUserInterrupt();	// check for user interrupt
      
	// interpolate the covar functions for the covariates
	table_lookup(&covariate_table,*time,cp,0);
      
	for (j = 0; j < nreps; j++, ft += nvars) { // loop over replicates
	
	  xp = &xs[nvars*((j%nrepx)+nrepx*k)];
	  pp = &ps[npars*(j%nrepp)];
	
	  (*ff)(ft,xp,pp,sidx,pidx,cidx,ncovars,cp,*time);
	
	}
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
