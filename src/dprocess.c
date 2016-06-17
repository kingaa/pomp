// dear emacs, please treat this as -*- C++ -*-

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include "pomp_internal.h"

SEXP do_dprocess (SEXP object, SEXP x, SEXP times, SEXP params, SEXP log, SEXP gnsi)
{
  int nprotect = 0;
  int *xdim, npars, nvars, nreps, nrepsx, ntimes;
  SEXP X, fn, fcall, rho;
  SEXP dimF;

  PROTECT(gnsi = duplicate(gnsi)); nprotect++;

  ntimes = length(times);
  if (ntimes < 2)
    error("dprocess error: length(times)<2: with no transitions, there is no work to do.");

  PROTECT(x = as_state_array(x)); nprotect++;
  xdim = INTEGER(GET_DIM(x)); 
  nvars = xdim[0]; nrepsx = xdim[1];
  if (ntimes != xdim[2])
    error("dprocess error: the length of 'times' and 3rd dimension of 'x' do not agree.");

  PROTECT(params = as_matrix(params)); nprotect++;
  xdim = INTEGER(GET_DIM(params)); 
  npars = xdim[0]; nreps = xdim[1];

  if (nrepsx > nreps) {         // more states than parameters
    if (nrepsx % nreps != 0) {
      error("dprocess error: the larger number of replicates is not a multiple of smaller.");
    } else {
      SEXP copy;
      double *src, *tgt;
      int dims[2];
      int j, k;
      dims[0] = npars; dims[1] = nrepsx;
      PROTECT(copy = duplicate(params)); nprotect++;
      PROTECT(params = makearray(2,dims)); nprotect++;
      setrownames(params,GET_ROWNAMES(GET_DIMNAMES(copy)),2);
      src = REAL(copy);
      tgt = REAL(params);
      for (j = 0; j < nrepsx; j++) {
        for (k = 0; k < npars; k++, tgt++) {
          *tgt = src[k+npars*(j%nreps)];
        }
      }
    }
    nreps = nrepsx;
  } else if (nrepsx < nreps) {  // more parameters than states
    if (nreps % nrepsx != 0) {
      error("dprocess error: the larger number of replicates is not a multiple of smaller.");
    } else {
      SEXP copy;
      double *src, *tgt;
      int dims[3];
      int i, j, k;
      dims[0] = nvars; dims[1] = nreps; dims[2] = ntimes;
      PROTECT(copy = duplicate(x)); nprotect++;
      PROTECT(x = makearray(3,dims)); nprotect++;
      setrownames(x,GET_ROWNAMES(GET_DIMNAMES(copy)),3);
      src = REAL(copy);
      tgt = REAL(x);
      for (i = 0; i < ntimes; i++) {
	for (j = 0; j < nreps; j++) {
	  for (k = 0; k < nvars; k++, tgt++) {
	    *tgt = src[k+nvars*((j%nrepsx)+nrepsx*i)];
	  }
	}
      }
    }
  }

  // extract the process function
  PROTECT(fn = GET_SLOT(object,install("dprocess"))); nprotect++;
  // construct the call
  PROTECT(fcall = VectorToPairList(GET_SLOT(object,install("userdata")))); nprotect++;
  PROTECT(fcall = LCONS(gnsi,fcall)); nprotect++;
  SET_TAG(fcall,install(".getnativesymbolinfo"));
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
    error("dprocess error: user 'dprocess' must return a rank-2 array.");
  }
  xdim = INTEGER(dimF);
  if ((xdim[0] != nreps) || (xdim[1] != ntimes-1)) {
    UNPROTECT(nprotect);
    error("dprocess error: user 'dprocess' must return a %d x %d array.",nreps,ntimes-1);
  }
  {
    const char *dimnms[2] = {"rep","time"};
    fixdimnames(X,dimnms,2);
  }
  UNPROTECT(nprotect);
  return X;
}
