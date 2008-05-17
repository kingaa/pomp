// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

SEXP do_dmeasure (SEXP object, SEXP y, SEXP x, SEXP times, SEXP params, SEXP log)
{
  int nprotect = 0;
  int *xdim, nvars, nreps, ntimes;
  SEXP X, fn, slotname, userdata, fcall, rho;
  SEXP dimP, dimX, dimY, dimF;
  ntimes = length(times);
  if (ntimes < 1)
    error("dmeasure error: no work to do");
  PROTECT(dimY = GET_DIM(y)); nprotect++;
  if ((isNull(dimY)) || (length(dimY)!=2))
    error("dmeasure error: 'y' must be a rank-2 array");
  PROTECT(dimX = GET_DIM(x)); nprotect++;
  if ((isNull(dimX)) || (length(dimX)!=3))
    error("dmeasure error: 'x' must be a rank-3 array");
  PROTECT(dimP = GET_DIM(params)); nprotect++;
  if ((isNull(dimP)) || (length(dimP)!=2))
    error("dmeasure error: 'params' must be a rank-2 array");
  xdim = INTEGER(dimX);
  nvars = xdim[0]; nreps = xdim[1];
  if (nreps != INTEGER(dimP)[1])
    error("dmeasure error: 2nd dimensions of 'params' and 'x' do not agree");
  if (ntimes != INTEGER(dimX)[2])
    error("rprocess error: length of 'times' and 3rd dimension of 'x' do not agree");
  if (ntimes != INTEGER(dimY)[1])
    error("rprocess error: length of 'times' and 2nd dimension of 'y' do not agree");
  PROTECT(slotname = NEW_CHARACTER(1)); nprotect++;
  // extract the process function
  SET_STRING_ELT(slotname,0,mkChar("dmeasure"));
  PROTECT(fn = GET_SLOT(object,slotname)); nprotect++;
  // extract the userda)ta
  SET_STRING_ELT(slotname,0,mkChar("userdata"));
  PROTECT(userdata = GET_SLOT(object,slotname)); nprotect++;
  // construct the call
  PROTECT(fcall = LCONS(fn,LCONS(y,LCONS(x,LCONS(times,LCONS(params,LCONS(log,VectorToPairList(userdata)))))))); nprotect++;
  PROTECT(rho = (CLOENV(fn))); nprotect++; // environment of the function
  PROTECT(X = eval(fcall,rho)); nprotect++; // do the call
  PROTECT(dimF = GET_DIM(X)); nprotect++;
  if ((isNull(dimF)) || (length(dimF) != 2))
    error("dmeasure error: user 'dmeasure' must return a rank-2 array");
  xdim = INTEGER(dimF);
  if ((xdim[0] != nreps) || (xdim[1] != ntimes))
    error("dmeasure error: user 'dmeasure' must return a %d x %d array",nreps,ntimes);
  UNPROTECT(nprotect);
  return X;
}

