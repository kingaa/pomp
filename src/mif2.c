// -*- C++ -*-

#include "pomp_internal.h"
#include <Rdefines.h>

SEXP randwalk_perturbation (SEXP params, SEXP rw_sd)
{
  int nprotect = 0;
  double *xp = 0, *rw, *xrw, *xs;
  SEXP Pnames, rwnames, pindex;
  int *dim, *pidx;
  int nrw = 0, npars, nreps;
  int j, k;

  // unpack parameter matrix
  PROTECT(params = duplicate(params)); nprotect++;
  xp = REAL(params);
  dim = INTEGER(GET_DIM(params)); npars = dim[0]; nreps = dim[1];
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;

  // names of parameters undergoing random walk
  PROTECT(rwnames = GET_NAMES(rw_sd)); nprotect++; 
  nrw = LENGTH(rwnames); rw = REAL(rw_sd);

  // indices of parameters undergoing random walk
  PROTECT(pindex = matchnames(Pnames,rwnames,"parameters")); nprotect++; 
  pidx = INTEGER(pindex);

  GetRNGstate();

  for (j = 0, xrw = rw; j < nrw; j++, pidx++, xrw++) {
    for (k = 0, xs = xp+(*pidx); k < nreps; k++, xs += npars) {
      *xs += rnorm(0,*xrw); 
    }
  }

  PutRNGstate();

  UNPROTECT(nprotect);
  return(params);
}
