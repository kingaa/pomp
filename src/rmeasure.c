// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

SEXP do_rmeasure (SEXP object, SEXP x, SEXP times, SEXP params)
{
  int nprotect = 0;
  int *xdim, nvars, nreps, ntimes, nobs;
  SEXP X, fn, slotname, userdata, fcall, rho;
  SEXP dimP, dimX, dimF, dimD;
  ntimes = length(times);
  if (ntimes < 1)
    error("rmeasure error: no work to do");
  PROTECT(dimX = GET_DIM(x)); nprotect++;
  if ((isNull(dimX)) || (length(dimX)!=3))
    error("rmeasure error: 'x' must be a rank-3 array");
  PROTECT(dimP = GET_DIM(params)); nprotect++;
  if ((isNull(dimP)) || (length(dimP)!=2))
    error("rmeasure error: 'params' must be a rank-2 array");
  xdim = INTEGER(dimX);
  nvars = xdim[0]; nreps = xdim[1];
  if (nreps != INTEGER(dimP)[1])
    error("rmeasure error: 2nd dimensions of 'params' and 'x' do not agree");
  if (ntimes != INTEGER(dimX)[2])
    error("rprocess error: length of 'times' and 3rd dimension of 'x' do not agree");
  PROTECT(slotname = NEW_CHARACTER(1)); nprotect++;
  // extract the process function
  SET_STRING_ELT(slotname,0,mkChar("rmeasure"));
  PROTECT(fn = GET_SLOT(object,slotname)); nprotect++;
  // extract the userda)ta
  SET_STRING_ELT(slotname,0,mkChar("userdata"));
  PROTECT(userdata = GET_SLOT(object,slotname)); nprotect++;
  // construct the call
  PROTECT(fcall = LCONS(fn,LCONS(x,LCONS(times,LCONS(params,VectorToPairList(userdata)))))); nprotect++;
  PROTECT(rho = (CLOENV(fn))); nprotect++; // environment of the function
  PROTECT(X = eval(fcall,rho)); nprotect++; // do the call
  PROTECT(dimF = GET_DIM(X)); nprotect++;
  if ((isNull(dimF)) || (length(dimF) != 3))
    error("rmeasure error: user 'rmeasure' must return a rank-3 array");
  SET_STRING_ELT(slotname,0,mkChar("data"));
  PROTECT(dimD = GET_DIM(GET_SLOT(object,slotname))); nprotect++;
  if ((isNull(dimD)) || (length(dimD) != 2)) // should not be necessary
    error("rmeasure error: slot 'data' must be a rank-2 array");
  nobs = INTEGER(dimD)[0];	// do we really need to check this?
  xdim = INTEGER(dimF);
  if ((xdim[0] != nobs) || (xdim[1] != nreps) || (xdim[2] != ntimes))
    error("rmeasure error: user 'rmeasure' must return a %d x %d x %d array",nobs,nreps,ntimes);
  if (isNull(GET_ROWNAMES(GET_DIMNAMES(X))))
    error("rmeasure error: user 'rmeasure' must return an array with rownames");
  UNPROTECT(nprotect);
  return X;
}

