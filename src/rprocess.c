// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

SEXP do_rprocess (SEXP object, SEXP xstart, SEXP times, SEXP params)
{
  int nprotect = 0;
  int *xdim, nvars, nreps, ntimes;
  SEXP X, fn, slotname, userdata, fcall, rho;
  SEXP dimXstart, dimP, dimX;
  ntimes = length(times);
  if (ntimes < 2)
    error("rprocess error: no transitions, no work to do");
  PROTECT(dimXstart = GET_DIM(xstart)); nprotect++;
  if ((isNull(dimXstart)) || (length(dimXstart)!=2))
    error("rprocess error: 'xstart' must be a rank-2 array");
  PROTECT(dimP = GET_DIM(params)); nprotect++;
  if ((isNull(dimP)) || (length(dimP)!=2))
    error("rprocess error: 'params' must be a rank-2 array");
  xdim = INTEGER(dimXstart);
  nvars = xdim[0]; nreps = xdim[1];
  if (nreps != INTEGER(dimP)[1])
    error("rprocess error: number of columns of 'params' and 'xstart' do not agree");
  PROTECT(slotname = NEW_CHARACTER(1)); nprotect++;
  // extract the process function
  SET_STRING_ELT(slotname,0,mkChar("rprocess"));
  PROTECT(fn = GET_SLOT(object,slotname)); nprotect++;
  // extract the userdata
  SET_STRING_ELT(slotname,0,mkChar("userdata"));
  PROTECT(userdata = GET_SLOT(object,slotname)); nprotect++;
  // construct the call
  PROTECT(fcall = LCONS(fn,LCONS(xstart,LCONS(times,LCONS(params,VectorToPairList(userdata)))))); nprotect++;
  PROTECT(rho = (CLOENV(fn))); nprotect++; // environment of the function
  PROTECT(X = eval(fcall,rho)); nprotect++; // do the call
  PROTECT(dimX = GET_DIM(X)); nprotect++;
  if ((isNull(dimX)) || (length(dimX) != 3))
    error("rprocess error: user 'rprocess' must return a rank-3 array");
  xdim = INTEGER(dimX);
  if ((xdim[0] != nvars) || (xdim[1] != nreps) || (xdim[2] != ntimes))
    error("rprocess error: user 'rprocess' must return a %d x %d x %d array",nvars,nreps,ntimes);
  if (isNull(GET_ROWNAMES(GET_DIMNAMES(X))))
    error("rprocess error: user 'rprocess' must return an array with rownames");
  UNPROTECT(nprotect);
  return X;
}
