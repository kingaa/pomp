// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include "pomp_internal.h"

SEXP do_partrans (SEXP object, SEXP params, SEXP dir, SEXP gnsi)
{
  int nprotect = 0;
  SEXP fn, fcall, rho, ans, nm;
  SEXP pdim, pvec;
  SEXP pompfun;
  SEXP tparams = R_NilValue;
  pompfunmode mode = undef;
  char direc;
  int qmat;
  int ndim[2], *dim, *idx;
  double *pp, *ps, *pt, *pa;
  int npar1, npar2, nreps;
  pomp_transform_fn *ff = NULL;
  int k;

  direc = *(INTEGER(dir));
  // extract the user-defined function
  switch (direc) {
  case -1: default:	// forward transformation
    PROTECT(pompfun = GET_SLOT(GET_SLOT(object,install("partrans")),
      install("from"))); nprotect++;
      PROTECT(fn = pomp_fun_handler(pompfun,gnsi,&mode)); nprotect++;
      break;
  case 1:			// inverse transformation
    PROTECT(pompfun = GET_SLOT(GET_SLOT(object,install("partrans")),
      install("to"))); nprotect++;
      PROTECT(fn = pomp_fun_handler(pompfun,gnsi,&mode)); nprotect++;
      break;
  break;
  }

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

  case Rfun: 			// use user-supplied R function

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
      const char *dimnm[2] = {"variable","rep"};
      ps = REAL(params);
      pp = REAL(pvec);

      memcpy(pp,ps,npar1*sizeof(double));

      PROTECT(ans = eval(fcall,rho)); nprotect++;

      PROTECT(nm = GET_NAMES(ans)); nprotect++;
      if (isNull(nm))
        errorcall(R_NilValue,"user transformation functions must return a named numeric vector");

      // set up matrix to hold the results
      npar2 = LENGTH(ans);
      ndim[0] = npar2; ndim[1] = nreps;
      PROTECT(tparams = makearray(2,ndim)); nprotect++;
      setrownames(tparams,nm,2);
      fixdimnames(tparams,dimnm,2);
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
        errorcall(R_NilValue,"user transformation functions must return a named numeric vector");

    }

    break;

  case native:			// use native routine

    *((void **) (&ff)) = R_ExternalPtrAddr(fn);

    if (qmat) {
      idx = INTEGER(PROTECT(name_index(GET_ROWNAMES(GET_DIMNAMES(params)),pompfun,"paramnames","parameters"))); nprotect++;
    } else {
      idx = INTEGER(PROTECT(name_index(PROTECT(GET_NAMES(params)),pompfun,"paramnames","parameters"))); nprotect+=2;
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

    UNPROTECT(nprotect);
    return params;

    break;

  }

  UNPROTECT(nprotect);
  return tparams;
}
