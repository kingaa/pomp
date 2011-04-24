// -*- C++ -*-

#include "pomp_internal.h"
#include <Rdefines.h>

static void nosort_resamp (int n, double *w, int *p, int offset) 
{
  int i, j;
  double du, u;

  for (j = 1; j < n; j++) w[j] += w[j-1];

  if (w[n-1] <= 0.0)
    error("non-positive sum of weights");

  GetRNGstate();
  du = w[n-1] / ((double) n);
  u = runif(-du,0);
  PutRNGstate();

  for (i = 0, j = 0; j < n; j++) {
    u += du;
    while (u > w[i]) i++;
    p[j] = i;
  }
  if (offset)			// add offset if needed
    for (j = 0; j < n; j++) p[j] += offset;

}

SEXP systematic_resampling (SEXP weights)
{
  int nprotect = 0;
  double u, du, *w;
  int i, j, *p, n;
  SEXP perm;

  n = LENGTH(weights);
  PROTECT(perm = NEW_INTEGER(n)); nprotect++;
  nosort_resamp(n,REAL(weights),INTEGER(perm),1);
  UNPROTECT(nprotect);
  return(perm);
}

SEXP pfilter_resample (SEXP weights, SEXP states, SEXP params) 
{
  int nprotect = 0;
  int xdim[2], nvars, npars, nreps;
  SEXP sdim, pdim, newstates, newparams, sample;
  SEXP retval, retvalnames;
  int i, j, n;
  double *ss, *st, *ps, *pt, *xx;

  PROTECT(sdim = GET_DIM(states)); nprotect++;
  nvars = INTEGER(sdim)[0]; nreps = INTEGER(sdim)[1];

  xdim[0] = nvars; xdim[1] = nreps;
  PROTECT(newstates = makearray(2,xdim)); nprotect++;
  setrownames(newstates,GET_ROWNAMES(GET_DIMNAMES(states)),2);
  ss = REAL(states);
  st = REAL(newstates);

  PROTECT(pdim = GET_DIM(params)); nprotect++;
  npars = INTEGER(pdim)[0]; 
  if (nreps != INTEGER(pdim)[1]) 
    error("'states' and 'params' do not agree in second dimension");

  xdim[0] = npars; xdim[1] = nreps;
  PROTECT(newparams = makearray(2,xdim)); nprotect++;
  setrownames(newparams,GET_ROWNAMES(GET_DIMNAMES(params)),2);
  ps = REAL(params);
  pt = REAL(newparams);

  {
    int sample[nreps];
    nosort_resamp(nreps,REAL(weights),sample,0);
    for (i = 0; i < nreps; i++) {
      for (j = 0, xx = ss+nvars*sample[i]; j < nvars; j++, st++, xx++) *st = *xx;
      for (j = 0, xx = ps+npars*sample[i]; j < npars; j++, pt++, xx++) *pt = *xx;
    }
  }

  PROTECT(retval = NEW_LIST(2)); nprotect++;
  SET_ELEMENT(retval,0,newstates);
  SET_ELEMENT(retval,1,newparams);
  PROTECT(retvalnames = NEW_CHARACTER(2)); nprotect++;
  SET_STRING_ELT(retvalnames,0,mkChar("states"));
  SET_STRING_ELT(retvalnames,1,mkChar("params"));
  SET_NAMES(retval,retvalnames);
  
  UNPROTECT(nprotect);
  return retval;
}
