// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include "pomp_internal.h"

SEXP do_init_state (SEXP object, SEXP params, SEXP t0)
{
  int nprotect = 0;
  SEXP x, x1, par, fcall, fn, rho, covar, tcovar, covars;
  int *dim, npar, nrep, nvar, index = 0;
  int xdim[2], j, k;
  double *p, *pp, *xp, *xpp;

  dim = INTEGER(GET_DIM(params));
  npar = dim[0]; nrep = dim[1]; 
  p = REAL(params);

  PROTECT(par = NEW_NUMERIC(npar)); nprotect++;
  SET_NAMES(par,GET_ROWNAMES(GET_DIMNAMES(params)));
  pp = REAL(par); 

  // extract the initializer function and its environment
  PROTECT(fn = GET_SLOT(object,install("initializer"))); nprotect++;
  PROTECT(rho = (CLOENV(fn))); nprotect++;

  // begin to construct the call
  PROTECT(fcall = VectorToPairList(GET_SLOT(object,install("userdata")))); nprotect++;
  PROTECT(fcall = LCONS(GET_SLOT(object,install("covarnames")),fcall)); nprotect++;
  SET_TAG(fcall,install("covarnames"));
  PROTECT(fcall = LCONS(GET_SLOT(object,install("paramnames")),fcall)); nprotect++;
  SET_TAG(fcall,install("paramnames"));
  PROTECT(fcall = LCONS(GET_SLOT(object,install("statenames")),fcall)); nprotect++;
  SET_TAG(fcall,install("statenames"));

  // extract covariates and interpolate
  PROTECT(tcovar = GET_SLOT(object,install("tcovar"))); nprotect++;
  if (LENGTH(tcovar) > 0) {	// do table lookkup
    PROTECT(covar = GET_SLOT(object,install("covar"))); nprotect++;
    PROTECT(covars = lookup_in_table(tcovar,covar,t0,&index)); nprotect++;
    PROTECT(fcall = LCONS(covars,fcall)); nprotect++;
    SET_TAG(fcall,install("covars"));
  }

  // finish constructing the call
  PROTECT(fcall = LCONS(t0,fcall)); nprotect++;
  SET_TAG(fcall,install("t0"));
  PROTECT(fcall = LCONS(par,fcall)); nprotect++;
  SET_TAG(fcall,install("params"));
  PROTECT(fcall = LCONS(fn,fcall)); nprotect++;

  for (k = 0; k < npar; k++) pp[k] = p[k];
  PROTECT(x1 = eval(fcall,rho)); nprotect++; // do the call

  if (!IS_NUMERIC(x1) || isNull(GET_NAMES(x1))) {
    UNPROTECT(nprotect);
    error("init.state error: user 'initializer' must return a named numeric vector");
  }

  nvar = LENGTH(x1);
  xp = REAL(x1);

  xdim[0] = nvar; xdim[1] = nrep;
  PROTECT(x = makearray(2,xdim)); nprotect++;
  setrownames(x,GET_NAMES(x1),2);
  xpp = REAL(x);

  for (k = 0; k < nvar; k++) xpp[k] = xp[k];

  for (j = 1; j < nrep; j++) {
    for (k = 0; k < npar; k++) pp[k] = p[j*npar+k];
    PROTECT(x1 = eval(fcall,rho)); nprotect++;
    xp = REAL(x1);
    for (k = 0; k < nvar; k++) xpp[j*nvar+k] = xp[k];
  }

  UNPROTECT(nprotect);
  return x;
}
