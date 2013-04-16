// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include "pomp_internal.h"

SEXP do_partrans (SEXP object, SEXP params, SEXP fun)
{
  int nprotect = 0;
  SEXP fn, fcall, rho, ans, nm;
  SEXP pdim, pvec;
  SEXP tparams = R_NilValue;
  int mode = -1;
  int qmat;
  int ndim[2], *dim, *idx;
  double *pp, *ps, *pt, *pa;
  int npar1, npar2, nreps;
  pomp_transform_fn *ff = NULL;
  int k;

  // extract the user-defined function
  PROTECT(fn = unpack_pomp_fun(fun,&mode)); nprotect++;
  // extract 'userdata' as pairlist
  PROTECT(fcall = VectorToPairList(GET_SLOT(object,install("userdata")))); nprotect++;

  PROTECT(pdim = GET_DIM(params)); nprotect++;
  if (isNull(pdim)) {		// a single vector
    npar1 = LENGTH(params); nreps = 1;
    qmat = 0;
  } else {			// a parameter matrix
    dim = INTEGER(pdim);
    npar1 = dim[0]; nreps = dim[1];
    qmat = 1;
  }

  switch (mode) {

  case 0: 			// use user-supplied R function

    // set up the function call
    if (qmat) {		// matrix case
      PROTECT(pvec = NEW_NUMERIC(npar1)); nprotect++;
      SET_NAMES(pvec,GET_ROWNAMES(GET_DIMNAMES(params)));
      PROTECT(fcall = LCONS(pvec,fcall)); nprotect++;
    } else {			// vector case
      PROTECT(fcall = LCONS(params,fcall)); nprotect++;
    }
    SET_TAG(fcall,install("params"));
    PROTECT(fcall = LCONS(fn,fcall)); nprotect++;

    // the function's environment
    PROTECT(rho = (CLOENV(fn))); nprotect++;

    if (qmat) {		// matrix case
      ps = REAL(params);
      pp = REAL(pvec);

      memcpy(pp,ps,npar1*sizeof(double));

      PROTECT(ans = eval(fcall,rho)); nprotect++;

      PROTECT(nm = GET_NAMES(ans)); nprotect++;
      if (isNull(nm))
	error("user transformation functions must return a named numeric vector");
      
      // set up matrix to hold the results
      npar2 = LENGTH(ans);
      ndim[0] = npar2; ndim[1] = nreps;
      PROTECT(tparams = makearray(2,ndim)); nprotect++;
      setrownames(tparams,nm,2);
      pt = REAL(tparams);

      pa = REAL(AS_NUMERIC(ans));
      memcpy(pt,pa,npar2*sizeof(double));

      ps += npar1;
      pt += npar2;
      for (k = 1; k < nreps; k++, ps += npar1, pt += npar2) {
	memcpy(pp,ps,npar1*sizeof(double));
	pa = REAL(AS_NUMERIC(eval(fcall,rho)));
	memcpy(pt,pa,npar2*sizeof(double));
      }
      
    } else {			// vector case
      
      PROTECT(tparams = eval(fcall,rho)); nprotect++;
      if (isNull(GET_NAMES(tparams)))
	error("user transformation functions must return a named numeric vector");
      
    }

    break;

  case 1:			// use native routine

    ff = (pomp_transform_fn *) R_ExternalPtrAddr(fn);
    
    if (qmat) {
      idx = INTEGER(PROTECT(name_index(GET_ROWNAMES(GET_DIMNAMES(params)),object,"paramnames"))); nprotect++;
    } else {
      idx = INTEGER(PROTECT(name_index(GET_NAMES(params),object,"paramnames"))); nprotect++;
    }

    set_pomp_userdata(fcall);

    PROTECT(tparams = duplicate(params)); nprotect++;

    for (k = 0, ps = REAL(params), pt = REAL(tparams); k < nreps; k++, ps += npar1, pt += npar1) {
      R_CheckUserInterrupt();
      (*ff)(pt,ps,idx);
    }

    unset_pomp_userdata();

    break;

  default:
    error("unrecognized 'mode' slot in 'partrans'");
  }

  UNPROTECT(nprotect);
  return tparams;
}
