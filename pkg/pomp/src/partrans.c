// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include "pomp_internal.h"

static SEXP pomp_R_transform_params (SEXP object, SEXP params, SEXP fun) {
  int nprotect = 0;
  SEXP Dim, nm;
  SEXP params_transformed, par, fcall, rho, ans;
  double *p = 0, *pp = 0, *pt = 0, *ps = 0;
  int nreps, npar1, npar2;
  int k;

  PROTECT(Dim = GET_DIM(params)); nprotect++;
  if (isNull(Dim)) {		// a single vector
    nreps = 1;
  } else {			// a parameter matrix
    int *dim = INTEGER(Dim);
    npar1 = dim[0]; nreps = dim[1];
  }

  if (nreps > 1) {		// matrix case
    PROTECT(par = NEW_NUMERIC(npar1)); nprotect++;
    SET_NAMES(par,GET_ROWNAMES(GET_DIMNAMES(params)));
  }    

  // set up the function call
  PROTECT(rho = (CLOENV(fun))); nprotect++;
  PROTECT(fcall = VectorToPairList(GET_SLOT(object,install("userdata")))); nprotect++;
  if (nreps > 1) {		// matrix case
    PROTECT(fcall = LCONS(par,fcall)); nprotect++;
  } else {			// vector case
    PROTECT(fcall = LCONS(params,fcall)); nprotect++;
  }
  SET_TAG(fcall,install("params"));
  PROTECT(fcall = LCONS(fun,fcall)); nprotect++;

  if (nreps > 1) {
    p = REAL(params);
    pp = REAL(par);

    memcpy(pp,p,npar1*sizeof(double));

    PROTECT(ans = eval(fcall,rho)); nprotect++;

    PROTECT(nm = GET_NAMES(ans)); nprotect++;
    if (isNull(nm))
      error("user transformation functions must return a named numeric vector");

    {
      int ndim[2];
      npar2 = LENGTH(ans);
      ndim[0] = npar2; ndim[1] = nreps;
      PROTECT(params_transformed = makearray(2,ndim)); nprotect++;
      setrownames(params_transformed,nm,2);
    }

    ps = REAL(AS_NUMERIC(ans));
    pt = REAL(params_transformed);
    memcpy(pt,ps,npar2*sizeof(double));

    p += npar1;
    pt += npar2;
    for (k = 1; k < nreps; k++, p += npar1, pt += npar2) {
      memcpy(pp,p,npar1*sizeof(double));
      //      PROTECT(ans = AS_NUMERIC(eval(fcall,rho))); 
      //      ps = REAL(ans);
      ps = REAL(AS_NUMERIC(eval(fcall,rho)));
      memcpy(pt,ps,npar2*sizeof(double));
      //      UNPROTECT(1);
    }

  } else {

    PROTECT(params_transformed = eval(fcall,rho)); nprotect++;
    if (isNull(GET_NAMES(params_transformed)))
      error("user transformation functions must return a named numeric vector");

  }

  UNPROTECT(nprotect);
  return params_transformed;
}

SEXP do_partrans (SEXP object, SEXP params, SEXP fun)
{
  int nprotect = 0;
  SEXP ptrans, fn, userdata;
  int use, drop;

  PROTECT(fn = VECTOR_ELT(fun,0)); nprotect++;
  use = INTEGER(VECTOR_ELT(fun,1))[0];

  switch (use) {
  case 0: 			// use user-supplied R function
    PROTECT(ptrans = pomp_R_transform_params(object,params,fn)); nprotect++;
    break;
  case 1:			// use native routine
    {
      pomp_transform_fn *ff = (pomp_transform_fn *) R_ExternalPtrAddr(fn);
      SEXP paramnames, pindex, Dim;
      int *idx, npar, nrep, qvec;
      double *ps, *pt;
      int k;

      PROTECT(Dim = GET_DIM(params)); nprotect++;
      if (isNull(Dim)) {	// a single vector
	npar = LENGTH(Dim); nrep = 1;
	qvec = 1;
      } else {			// a parameter matrix
	int *dim;
	dim = INTEGER(Dim);
	npar = dim[0]; nrep = dim[1];
	qvec = 0;
      }

      PROTECT(paramnames = GET_SLOT(object,install("paramnames"))); nprotect++;
      if (LENGTH(paramnames) > 0) {
	if (qvec) {
	  PROTECT(pindex = matchnames(GET_NAMES(params),paramnames)); nprotect++;
	} else {
	  PROTECT(pindex = matchnames(GET_ROWNAMES(GET_DIMNAMES(params)),paramnames)); nprotect++;
	}
	idx = INTEGER(pindex);
      } else {
	idx = 0;
      }

      PROTECT(userdata = VectorToPairList(GET_SLOT(object,install("userdata")))); nprotect++;
      set_pomp_userdata(userdata);

      PROTECT(ptrans = duplicate(params)); nprotect++;

      for (k = 0, ps = REAL(params), pt = REAL(ptrans); k < nrep; k++, ps += npar, pt += npar) {
	R_CheckUserInterrupt();
	(*ff)(pt,ps,idx);
      }

      unset_pomp_userdata();

    }
    break;
  default:
    error("unrecognized 'use' slot in 'partrans'");
  }

  UNPROTECT(nprotect);
  return ptrans;
}
