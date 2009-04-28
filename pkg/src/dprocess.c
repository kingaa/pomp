// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

SEXP do_dprocess (SEXP object, SEXP x, SEXP times, SEXP params, SEXP log)
{
  int nprotect = 0;
  int *xdim, nvars, nreps, ntimes;
  SEXP X, fn, fcall, rho;
  SEXP dimP, dimX, dimF;
  ntimes = length(times);
  if (ntimes < 2) {
    UNPROTECT(nprotect);
    error("dprocess error: no transitions, no work to do");
  }
  PROTECT(dimX = GET_DIM(x)); nprotect++;
  if ((isNull(dimX)) || (length(dimX)!=3)) {
    UNPROTECT(nprotect);
    error("dprocess error: 'x' must be a rank-3 array");
  }
  PROTECT(dimP = GET_DIM(params)); nprotect++;
  if ((isNull(dimP)) || (length(dimP)!=2)) {
    UNPROTECT(nprotect);
    error("dprocess error: 'params' must be a rank-2 array");
  }
  xdim = INTEGER(dimX);
  nvars = xdim[0]; nreps = xdim[1];
  if (nreps != INTEGER(dimP)[1]) {
    UNPROTECT(nprotect);
    error("dprocess error: number of columns of 'params' and 'x' do not agree");
  }
  if (ntimes != INTEGER(dimX)[2]) {
    UNPROTECT(nprotect);
    error("rprocess error: length of 'times' and 3rd dimension of 'x' do not agree");
  }
  // extract the process function
  PROTECT(fn = GET_SLOT(object,install("dprocess"))); nprotect++;
  // construct the call
  PROTECT(fcall = VectorToPairList(GET_SLOT(object,install("userdata")))); nprotect++;
  PROTECT(fcall = LCONS(GET_SLOT(object,install("PACKAGE")),fcall)); nprotect++;
  SET_TAG(fcall,install("PACKAGE"));
  PROTECT(fcall = LCONS(GET_SLOT(object,install("covarnames")),fcall)); nprotect++;
  SET_TAG(fcall,install("covarnames"));
  PROTECT(fcall = LCONS(GET_SLOT(object,install("paramnames")),fcall)); nprotect++;
  SET_TAG(fcall,install("paramnames"));
  PROTECT(fcall = LCONS(GET_SLOT(object,install("statenames")),fcall)); nprotect++;
  SET_TAG(fcall,install("statenames"));
  PROTECT(fcall = LCONS(GET_SLOT(object,install("covar")),fcall)); nprotect++;
  SET_TAG(fcall,install("covar"));
  PROTECT(fcall = LCONS(GET_SLOT(object,install("tcovar")),fcall)); nprotect++;
  SET_TAG(fcall,install("tcovar"));
  PROTECT(fcall = LCONS(log,fcall)); nprotect++;
  SET_TAG(fcall,install("log"));
  PROTECT(fcall = LCONS(params,fcall)); nprotect++;
  SET_TAG(fcall,install("params"));
  PROTECT(fcall = LCONS(AS_NUMERIC(times),fcall)); nprotect++;
  SET_TAG(fcall,install("times"));
  PROTECT(fcall = LCONS(x,fcall)); nprotect++;
  SET_TAG(fcall,install("x"));
  PROTECT(fcall = LCONS(fn,fcall)); nprotect++;
  PROTECT(rho = (CLOENV(fn))); nprotect++; // environment of the function
  PROTECT(X = eval(fcall,rho)); nprotect++; // do the call
  PROTECT(dimF = GET_DIM(X)); nprotect++;
  if ((isNull(dimF)) || (length(dimF) != 2)) {
    UNPROTECT(nprotect);
    error("dprocess error: user 'dprocess' must return a rank-2 array");
  }
  xdim = INTEGER(dimF);
  if ((xdim[0] != nreps) || (xdim[1] != ntimes-1)) {
    UNPROTECT(nprotect);
    error("dprocess error: user 'dprocess' must return a %d x %d array",nreps,ntimes-1);
  }
  UNPROTECT(nprotect);
  return X;
}

