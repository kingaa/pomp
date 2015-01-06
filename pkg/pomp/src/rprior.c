// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "pomp_internal.h"

//_pomp_default_rprior (double *p, int *parindex) {
//  return;
//}

SEXP do_rprior (SEXP object, SEXP params, SEXP gnsi)
{
  int nprotect = 0;
  int mode = -1;
  int npars, nreps;
  SEXP Pnames, P, fn, fcall;
  SEXP pompfun;
  int *dim;
  const char *dimnms[2] = {"variable","rep"};

  PROTECT(params = as_matrix(params)); nprotect++;
  dim = INTEGER(GET_DIM(params));
  npars = dim[0]; nreps = dim[1]; 

  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;
    
  // extract the user-defined function
  PROTECT(pompfun = GET_SLOT(object,install("rprior"))); nprotect++;
  PROTECT(fn = pomp_fun_handler(pompfun,gnsi,&mode)); nprotect++;

  // extract 'userdata' as pairlist
  PROTECT(fcall = VectorToPairList(GET_SLOT(object,install("userdata")))); nprotect++;

  // first do setup
  switch (mode) {
  case 0:			// use R function

    {
      SEXP pvec, rho, ans, nm;
      int first = 1;
      int use_names = 0;
      double *pa, *pp, *ps, *pt;
      int *posn;
      int i, j;

      // to store results
      PROTECT(P = makearray(2,dim)); nprotect++;
      setrownames(P,Pnames,2);
      
      // temporary storage
      PROTECT(pvec = NEW_NUMERIC(npars)); nprotect++;
      SET_NAMES(pvec,Pnames);
      
      // set up the function call
      PROTECT(fcall = LCONS(pvec,fcall)); nprotect++;
      SET_TAG(fcall,install("params"));
      PROTECT(fcall = LCONS(fn,fcall)); nprotect++;
      
      // get the function's environment
      PROTECT(rho = (CLOENV(fn))); nprotect++;
      
      pp = REAL(pvec);

      for (j = 0, ps = REAL(params), pt = REAL(P); j < nreps; j++, ps += npars, pt += npars) {

	for (i = 0; i < npars; i++) pp[i] = ps[i];

	if (first) {
	  // evaluate the call
	  PROTECT(ans = eval(fcall,rho)); nprotect++;
	  if (LENGTH(ans) != npars) {
	    error("user 'rprior' returns a vector of %d parameters but %d are expected",
		  LENGTH(ans),npars);
	  }
	  
	  // get name information to fix potential alignment problems
	  PROTECT(nm = GET_NAMES(ans)); nprotect++;
	  use_names = !isNull(nm);
	  if (use_names) {   // match names against names from params slot
	    posn = INTEGER(PROTECT(matchnames(Pnames,nm,"parameters"))); nprotect++;
	  } else {
	    posn = 0;
	  }
	  
	  pa = REAL(AS_NUMERIC(ans));
	  
	  first = 0;
	
	} else {
	  
	  pa = REAL(AS_NUMERIC(eval(fcall,rho)));

	}

	if (use_names) {
	  for (i = 0; i < npars; i++) pt[posn[i]] = pa[i];
	} else {
	  for (i = 0; i < npars; i++) pt[i] = pa[i];
	}
      }
    }

    break;

  case 1:			// use native routine

    {
      double *pp, *ps;
      int *pidx = 0;
      pomp_rprior *ff = NULL;
      int j;

      PROTECT(P = duplicate(params)); nprotect++;

      // construct state, parameter, covariate, observable indices
      pidx = INTEGER(PROTECT(name_index(Pnames,pompfun,"paramnames"))); nprotect++;
      
      // address of native routine
      ff = (pomp_rprior *) R_ExternalPtrAddr(fn);

      R_CheckUserInterrupt();	// check for user interrupt

      set_pomp_userdata(fcall);
      GetRNGstate();

      // loop over replicates
      for (j = 0, ps = REAL(P); j < nreps; j++, ps += npars)
	(*ff)(ps,pidx);

      PutRNGstate();
      unset_pomp_userdata();
    }
    
    break;

  default:

    error("unrecognized 'mode' slot in 'rprior'");
    break;

  }

  UNPROTECT(nprotect);
  fixdimnames(P,dimnms,2);
  return P;
}
