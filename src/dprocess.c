// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

SEXP do_dprocess (SEXP object, SEXP x, SEXP times, SEXP params, SEXP log)
{
  int nprotect = 0;
  int *xdim, nvars, nreps, ntimes;
  SEXP X, fn, slotname, userdata, fcall, rho;
  SEXP dimP, dimX, dimF;
  ntimes = length(times);
  if (ntimes < 2)
    error("dprocess error: no transitions, no work to do");
  PROTECT(dimX = GET_DIM(x)); nprotect++;
  if ((isNull(dimX)) || (length(dimX)!=3))
    error("dprocess error: 'x' must be a rank-3 array");
  PROTECT(dimP = GET_DIM(params)); nprotect++;
  if ((isNull(dimP)) || (length(dimP)!=2))
    error("dprocess error: 'params' must be a rank-2 array");
  xdim = INTEGER(dimX);
  nvars = xdim[0]; nreps = xdim[1];
  if (nreps != INTEGER(dimP)[1])
    error("dprocess error: number of columns of 'params' and 'x' do not agree");
  if (ntimes != INTEGER(dimX)[2])
    error("rprocess error: length of 'times' and 3rd dimension of 'x' do not agree");
  PROTECT(slotname = NEW_CHARACTER(1)); nprotect++;
  // extract the process function
  SET_STRING_ELT(slotname,0,mkChar("dprocess"));
  PROTECT(fn = GET_SLOT(object,slotname)); nprotect++;
  // extract the userda)ta
  SET_STRING_ELT(slotname,0,mkChar("userdata"));
  PROTECT(userdata = GET_SLOT(object,slotname)); nprotect++;
  // construct the call
  PROTECT(fcall = LCONS(fn,LCONS(x,LCONS(times,LCONS(params,LCONS(log,VectorToPairList(userdata))))))); nprotect++;
  PROTECT(rho = (CLOENV(fn))); nprotect++; // environment of the function
  PROTECT(X = eval(fcall,rho)); nprotect++; // do the call
  PROTECT(dimF = GET_DIM(X)); nprotect++;
  if ((isNull(dimF)) || (length(dimF) != 2))
    error("dprocess error: user 'dprocess' must return a rank-2 array");
  xdim = INTEGER(dimF);
  if ((xdim[0] != nreps) || (xdim[1] != ntimes-1))
    error("dprocess error: user 'dprocess' must return a %d x %d array",nreps,ntimes-1);
  UNPROTECT(nprotect);
  return X;
}

